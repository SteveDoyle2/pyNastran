# pylint: disable=R0902,R0904,R0914,C0111
"""
All rigid elements are defined in this file.  This includes:

 * RBAR
 * RBAR1
 * RBE1
 * RBE2
 * RBE3

All rigid elements are RigidElement and Element objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six import string_types, integer_types
from six.moves import zip, range
from itertools import count

from pyNastran.bdf.field_writer_8 import set_blank_if_default, print_card_8
from pyNastran.bdf.cards.baseCard import Element
from pyNastran.bdf.bdfInterface.assign_type import (integer,
    integer_or_double, integer_double_or_blank, integer_or_blank,
    double_or_blank, integer_double_or_string, components, components_or_blank,
    blank, string)
from pyNastran.bdf.field_writer_16 import print_card_16
# from pyNastran.bdf.cards.utils import build_table_lines


class RigidElement(Element):
    def cross_reference(self, model):
        pass


class RBAR(RigidElement):
    type = 'RBAR'

    def __init__(self, card=None, data=None, comment=''):
        """
        +------+-----+----+----+--------+-----+-----+-----+-------+
        | RBAR | EID | GA | GB | CNA    | CNB | CMA | CMB | ALPHA |
        +------+-----+----+----+--------+-----+-----+-----+-------+
        | RBAR | 5   | 1  |  2 | 123456 |     |     |     | 6.5-6 |
        +------+-----+----+----+--------+-----+-----+-----+-------+
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.ga = integer(card, 2, 'ga')
            self.gb = integer(card, 3, 'gb')
            self.cna = components_or_blank(card, 4, 'cna')
            self.cnb = components_or_blank(card, 5, 'cnb')
            self.cma = components_or_blank(card, 6, 'cma')
            self.cmb = components_or_blank(card, 7, 'cmb')
            self.alpha = double_or_blank(card, 8, 'alpha', 0.0)
            assert len(card) <= 9, 'len(RBAR card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.ga = data[1]
            self.gb = data[2]
            self.cna = data[3]
            self.cnb = data[4]
            self.cma = data[5]
            self.cmb = data[6]
            self.alpha = data[7]

    # def convert_to_MPC(self, mpcID):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError()
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpcID]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    #def write_code_aster(self):
        #msg = ''
        #msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBAR\n"
        #msg += "        MODELE=MODELE,\n"  # rigid element
        #msg += "        \n"
        #return msg

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, integer_types) or self.gb is None:
            return self.gb
        return self.gb.nid

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.Ga(), msg=msg)
        self.gb = model.Node(self.Gb(), msg=msg)

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        return [self.Ga()]

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return [self.Gb()]

    def raw_fields(self):
        list_fields = ['RBAR', self.eid, self.Ga(), self.Gb(), self.cna,
                       self.cnb, self.cma, self.cmb, self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR', self.eid, self.Ga(), self.Gb(), self.cna, self.cnb,
                       self.cma, self.cmb, alpha]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBAR1(RigidElement):
    type = 'RBAR1'

    def __init__(self, card=None, data=None, comment=''):
        """
        +-------+-----+----+----+-----+-------+
        | RBAR1 | EID | GA | GB | CB  | ALPHA |
        +-------+-----+----+----+-----+-------+
        | RBAR1 | 5   |  1 |  2 | 123 | 6.5-6 |
        +-------+-----+----+----+-----+-------+
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            self.eid = integer(card, 1, 'eid')
            self.ga = integer(card, 2, 'ga')
            self.gb = integer(card, 3, 'gb')
            self.cb = components_or_blank(card, 4, 'cb')
            self.alpha = double_or_blank(card, 5, 'alpha', 0.0)
            assert len(card) <= 6, 'len(RBAR1 card) = %i' % len(card)
        else:
            self.eid = data[0]
            self.ga = data[1]
            self.gb = data[2]
            self.cb = data[3]
            self.alpha = data[4]

    def Ga(self):
        if isinstance(self.ga, integer_types):
            return self.ga
        return self.ga.nid

    def Gb(self):
        if isinstance(self.gb, integer_types) or self.gb is None:
            return self.gb
        return self.gb.nid

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.ga = model.Node(self.Ga(), msg=msg)
        self.gb = model.Node(self.Gb(), msg=msg)

    def raw_fields(self):
        list_fields = ['RBAR1', self.eid, self.Ga(), self.Gb(), self.cb, self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.0)
        list_fields = ['RBAR1', self.eid, self.Ga(), self.Gb(), self.cb, alpha]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBE1(RigidElement):  # maybe not done, needs testing
    type = 'RBE1'

    def __init__(self, card=None, data=None, comment=''):
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment

        self.eid = integer(card, 1, 'eid')
        self.Gni = []
        self.Cni = []

        iUm = card.index('UM')
        if iUm > 0:
            assert string(card, iUm, 'UM') == 'UM'

        #assert isinstance(card[-1], string_types), 'card[-1]=%r type=%s' %(card[-1], type(card[-1]))
        alpha_last = integer_double_or_string(card, -1, 'alpha_last')
        if isinstance(alpha_last, float):
            self.alpha = alpha_last
            card.pop()  # remove the last field so len(card) will not include alpha
        else:
            self.alpha = 0.

        # loop till UM, no field9,field10
        n = 1
        i = 0
        offset = 2
        while offset + i < iUm - 1:
            #print('field(%s) = %s' % (offset + i, card.field(offset + i)))
            gni = integer_or_blank(card, offset + i, 'gn%i' % n)
            cni = components_or_blank(card, offset + i + 1, 'cn%i' % n)

            if gni:
                #print("gni=%s cni=%s" % (gni ,cni))
                self.Gni.append(gni)
                self.Cni.append(cni)
                n += 1
            else:
                assert cni is None
            i += 2

        self.Gmi = []
        self.Cmi = []

        # loop till alpha, no field9,field10
        n = 1
        offset = iUm + 1
        i = 0

        # dont grab alpha
        while offset + i < len(card):
            gmi = integer_or_blank(card, offset + i, 'gm%i' % n)
            cmi = components_or_blank(card, offset + i + 1, 'cm%i' % n)
            if gmi:
                self.Gmi.append(gmi)
                self.Cmi.append(cmi)
                n += 1
            else:
                assert cmi is None
            i += 2

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.Gni = model.Nodes(self.Gni, allowEmptyNodes=True, msg=msg)
        self.Gmi = model.Nodes(self.Gmi, allowEmptyNodes=True, msg=msg)

    @property
    def Gni_node_ids(self):
        if len(self.Gni) == 0:
            return []
        return self._nodeIDs(nodes=self.Gni, allowEmptyNodes=True)

    @property
    def Gmi_node_ids(self):
        if len(self.Gmi) == 0:
            return []
        return self._nodeIDs(nodes=self.Gmi, allowEmptyNodes=True)

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        # checked
        return self.Gni_node_ids

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        # checked
        nodes = self.Gmi_node_ids
        return nodes

    def raw_fields(self):
        list_fields = [self.type, self.eid]
        for (i, gn, cn) in zip(count(), self.Gni_node_ids, self.Cni):
            list_fields += [gn, cn]
            if i > 0 and i % 3 == 0:
                list_fields += [None]

        nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM/ALPHA onto next line
        if nspaces < 8:
            list_fields += [None] * nspaces

        # overly complicated loop to print the UM section
        list_fields += ['UM']
        j = 1
        for (i, gm, cm) in zip(count(), self.Gmi_node_ids, self.Cmi):
            list_fields += [gm, cm]
            if i > 0 and j % 3 == 0:
                list_fields += [None, None]
                j -= 3
            j += 1

        if self.alpha > 0.:  # handles default alpha value
            nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
            if nspaces == 1:
                list_fields += [None, None]
            list_fields += [self.alpha]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class RBE2(RigidElement):
    type = 'RBE2'
    _field_map = {1: 'eid', 2:'gn', 3:'cm'}

    def _update_field_helper(self, n, value):
        """
        Updates complicated parameters on the GRID card

        Parameters
        ----------
        self : RBE2()
            the RBE2 object pointer
        n : int
            the field number to update
        value : int/float/str
            the value for the appropriate field
        """
        if n > 3 and n <= 3 + len(self.Gmi):
            self.Gmi[n - 4] = value
        elif n == 4 + len(self.Gmi):
            self.alpha = value
        else:
            raise KeyError('Field %r is an invalid %s entry.' % (n, self.type))
        return value

    def __init__(self, card=None, data=None, comment=''):
        """
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |   1   |  2  |  3  |  4  |  5   |   6   |  7  |  8  |  9  |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |  RBE2 | EID | GN  | CM  | GM1  | GM2   | GM3 | GM4 | GM5 |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        |       | GM6 | GM7 | GM8 | etc. | ALPHA |     |     |     |
        +-------+-----+-----+-----+------+-------+-----+-----+-----+
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        if card:
            #: Element identification number
            self.eid = integer(card, 1, 'eid')

            #: Identification number of grid point to which all six independent
            #: degrees-of-freedom for the element are assigned. (Integer > 0)
            self.gn = integer(card, 2, 'gn')

            #: Component numbers of the dependent degrees-of-freedom in the
            #: global coordinate system at grid points GMi. (Integers 1 through
            #: 6 with no embedded blanks.)
            self.cm = components_or_blank(card, 3, 'cm')

            alpha = integer_or_double(card, len(card) - 1, 'alpha')
            if isinstance(alpha, float):
                #: Grid point identification numbers at which dependent
                #: degrees-of-freedom are assigned. (Integer > 0)
                self.alpha = alpha

                # the last field is not part of Gmi
                n = 1
            else:
                # the last field is part of Gmi
                n = 0
                self.alpha = 0.0

            j = 4
            self.Gmi = []
            for i in range(len(card) - 4 - n):
                gmi = integer(card, j + i, 'Gm%i' % (i + 1))
                self.Gmi.append(gmi)
        else:
            self.eid = data[0]
            self.gn = data[1]
            self.cm = data[2]
            self.Gmi = data[3]
            self.alpha = data[4]
            #print("eid=%s gn=%s cm=%s Gmi=%s alpha=%s"
                  #% (self.eid, self.gn, self.cm, self.Gmi, self.alpha))
            #raise NotImplementedError('RBE2 data...')

        assert self.gn is not None, 'gn=%s' % self.gn
        assert self.cm is not None, 'cm=%s' % self.cm
        self.gn = self.gn
        self.cm = str(self.cm)

    def convert_to_MPC(self, mpc_id):
        """
        .. math:: -A_i u_i + A_j u_j = 0

        where :math:`u_i` are the base DOFs (max=6)

         +------+------+----+----+-----+----+----+----+
         | MPC  | sid  | g1 | c1 | a1  | g2 | c2 | a2 |
         +------+------+----+----+-----+----+----+----+
         | RBE2 | eid  | gn | cm | g1  | g2 | g3 | g4 |
         +------+------+----+----+-----+----+----+----+
        """
        nCM = len(self.cm)
        Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1

        card = ['MPC', mpc_id]
        for cm in self.cm:
            # the minus sign is applied to the base node
            card += [self.gn, cm, -Ai]

        for gm in self.Gmi:
            for cm in self.cm:
                card += [gm, cm, Ai]
        return card

    #def convert_to_RBE3(self):
        #raise NotImplementedError()
        #eid = self.eid
        #ref_node = self.gn
        #dof = self.cm
        #wf = 1.0
        #sDof = 123  # this is probably wrong...
        #boundary_nodes = self.Gmi

        ## this is to get the farthest nodes for the UM card
        #boundary_nodes.sort()
        #rbe3_nodes = boundary_nodes

        #rbe3 = ['RBE3', eid, ref_node, dof, wf, sDof] + rbe3_nodes
        #return rbe3

    def write_code_aster(self):
        """
        Converts to a LIAISON SOLIDE for dofs 123456.
        For other dof combinations, general MPC equations are written
        """
        msg = ''
        msg += "BLOCAGE=AFFE_CHAR_MECA(  # RBE2 ID=%s\n" % (self.eid)
        msg += "        MODELE=MODELE,\n"  # rigid element
        if self.cm == 123456:
            msg += "        LIASON_SOLIDE=(\n"
            msg += "        _F(NOEUD=\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," % (nid)
            msg = msg[:-1]
            msg += '\n'
        else:
            msg += "        _F(NOEUD=  # doesnt handle coordinate systems\n"
            msg += "           "
            for nid in self.Gmi:
                msg += "'N%i'," % (nid)
            msg = msg[:-1]
            msg += '\n'
        return msg

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        self.Gmi = model.Nodes(self.Gmi_node_ids, allowEmptyNodes=True, msg=msg)
        self.gn = model.Node(self.Gn(), msg=msg)

    def Gn(self):
        if isinstance(self.gn, integer_types) or self.gn is None:
            return self.gn
        return self.gn.nid

    @property
    def Gmi_node_ids(self):
        if len(self.Gmi) == 0:
            return []
        return self._nodeIDs(nodes=self.Gmi, allowEmptyNodes=True)

    @property
    def independent_nodes(self):
        """gets the independent node ids"""
        nodes = [self.Gn()]
        return nodes

    @property
    def dependent_nodes(self):
        """gets the dependent node ids"""
        return self.Gmi_node_ids

    def raw_fields(self):
        list_fields = ['RBE2', self.eid, self.Gn(), self.cm] + self.Gmi_node_ids + [self.alpha]
        return list_fields

    def repr_fields(self):
        alpha = set_blank_if_default(self.alpha, 0.)
        list_fields = ['RBE2', self.eid, self.Gn(), self.cm] + self.Gmi_node_ids + [alpha]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)


class RBE3(RigidElement):
    """
    .. todo:: not done, needs testing badly
    """
    type = 'RBE3'

    def __init__(self, card=None, data=None, comment=''):
        """
        eid
        refgrid
        refc
        WtCG_groups = [wt, ci, Gij]
        Gmi
        Cmi
        alpha
        """
        RigidElement.__init__(self, card, data)
        if comment:
            self._comment = comment
        self.eid = integer(card, 1, 'eid')
        blank(card, 2, 'blank')
        self.refgrid = integer(card, 3, 'refgrid')
        self.refc = components_or_blank(card, 4, 'refc')

        fields = [field.upper() if isinstance(field, string_types) else field for field in card[5:]]
        ioffset = 5
        iwt_max = len(fields) + ioffset
        try:
            ialpha = fields.index('ALPHA') + ioffset
            iwt_max = ialpha  # the index to start parsing UM
            ium_stop = ialpha  # the index to stop  parsing UM
        except ValueError:
            ialpha = None
            ium_stop = iwt_max

        try:
            ium = fields.index('UM') + ioffset
            iwt_max = ium
        except ValueError:
            ium = None

        self.WtCG_groups = []

        i = ioffset
        n = 1
        while i < iwt_max:
            Gij = []
            wtname = 'wt' + str(n)
            wt = double_or_blank(card, i, wtname)
            if wt is not None:
                cname = 'c'+str(n)
                compi = components_or_blank(card, i + 1, cname)

                #print("%s=%s %s=%s" % (wtname, wt, cname, compi))
                i += 2
                gij = 0

                j = 0
                while isinstance(gij, int) and i < iwt_max:
                    j += 1
                    gij_name = 'g%s,%s' % (n, j)
                    gij = integer_double_or_blank(card, i, gij_name)
                    if isinstance(gij, float):
                        break
                    #print("%s = %s" % (gij_name, gij))
                    if gij is not None:
                        Gij.append(gij)
                    i += 1
                assert compi is not None
                assert len(Gij) > 0, Gij
                assert Gij[0] is not None, Gij
                wtCG_group = [wt, compi, Gij]
                self.WtCG_groups.append(wtCG_group)
                #print('----finished a group=%r----' % wtCG_group)
            else:
                i += 1

        self.Gmi = []
        self.Cmi = []
        if ium:
            #print('UM = %s' % card.field(ium))  # UM
            i = ium + 1
            n = 1
            #print("i=%s iUmStop=%s" % (i, iUmStop))
            for j in range(i, ium_stop, 2):

                gm_name = 'gm' + str(n)
                cm_name = 'cm' + str(n)
                gmi = integer_or_blank(card, j, gm_name)
                if gmi is not None:
                    cmi = components(card, j + 1, cm_name)
                    #print "gmi=%s cmi=%s" % (gmi, cmi)
                    self.Gmi.append(gmi)
                    self.Cmi.append(cmi)

        if ialpha:
            self.alpha = double_or_blank(card, ialpha + 1, 'alpha')
        else:
            #: thermal expansion coefficient
            self.alpha = 0.0

    # def convert_to_MPC(self, mpcID):
    #     """
    #     -Ai*ui + Aj*uj = 0
    #     where ui are the base DOFs (max=6)
    #     mpc sid   g1 c1 a1  g2 c2 a2
    #     rbe2 eid  gn cm g1  g2 g3 g4
    #     """
    #     raise NotImplementedError('this is the code for an RBE2...not RBE3')
    #     #i = 0
    #     nCM = len(self.cm)
    #     Ai = nCM * len(self.Gmi) / len(self.gn)  # where nGN=1
    #
    #     card = ['MPC', mpcID]
    #     for cm in self.cm:  # the minus sign is applied to the base node
    #         card += [self.gn, cm, -Ai]
    #
    #     for gm in self.Gmi:
    #         for cm in self.cm:
    #             card += [gm, cm, Ai]
    #     return card

    @property
    def ref_grid_id(self):
        if isinstance(self.refgrid, int) or self.refgrid is None:
            return self.refgrid
        return self.refgrid.nid

    @property
    def Gmi_node_ids(self):
        if len(self.Gmi) == 0:
            return []
        return self._nodeIDs(nodes=self.Gmi, allowEmptyNodes=True)

    def cross_reference(self, model):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        assert self.Gmi is not None
        self.Gmi = model.Nodes(self.Gmi, allowEmptyNodes=True, msg=msg)
        assert self.Gmi is not None
        self.refgrid = model.Node(self.ref_grid_id, msg=msg)
        for i, (wt, ci, Gij) in enumerate(self.WtCG_groups):
            self.WtCG_groups[i][2] = model.Nodes(Gij, allowEmptyNodes=True, msg=msg)

    def safe_cross_reference(self, model, debug=True):
        msg = ' which is required by %s eid=%s' % (self.type, self.eid)
        assert self.Gmi is not None
        self.Gmi = model.Nodes(self.Gmi, allowEmptyNodes=True, msg=msg)
        assert self.Gmi is not None
        self.refgrid = model.Node(self.ref_grid_id, msg=msg)
        for i, (wt, ci, Gij) in enumerate(self.WtCG_groups):
            self.WtCG_groups[i][2] = model.Nodes(Gij, allowEmptyNodes=True, msg=msg)

    @property
    def independent_nodes(self):
        """
        gets the independent node ids
        TODO: not checked
        """
        nodes = []
        for (wt, compi, Gij) in self.WtCG_groups:
            Giji = self._nodeIDs(nodes=Gij, allowEmptyNodes=True)
            nodes += Giji
        return nodes

    @property
    def dependent_nodes(self):
        """
        gets the dependent node ids
        TODO: not checked
        """
        nodes = [self.ref_grid_id]
        nodes += self.Gmi_node_ids
        return nodes

    def raw_fields(self):
        list_fields = ['RBE3', self.eid, None, self.ref_grid_id, self.refc]
        for (wt, ci, Gij) in self.WtCG_groups:
            Giji = self._nodeIDs(nodes=Gij, allowEmptyNodes=True)
            list_fields += [wt, ci] + Giji
        nspaces = 8 - (len(list_fields) - 1) % 8  # puts UM onto next line

        if nspaces < 8:
            list_fields += [None] * nspaces

        if self.Gmi:
            list_fields += ['UM']
            for (gmi, cmi) in zip(self.Gmi_node_ids, self.Cmi):
                list_fields += [gmi, cmi]

        nspaces = 8 - (len(list_fields) - 1) % 8  # puts ALPHA onto next line
        if nspaces < 8:
            list_fields += [None] * nspaces

        if self.alpha > 0.:  # handles the default value
            list_fields += ['ALPHA', self.alpha]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_8(card)
