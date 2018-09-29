# pylint: disable=C0103,R0902,R0904,R0914
"""
All dynamic control cards are defined in this file.  This includes:

 * FREQ
 * FREQ1
 * FREQ2 (not implemented)
 * FREQ3
 * FREQ4
 * FREQ5 (not implemented)
 * NLPCI
 * NLPARM
 * TSTEP
 * TSTEPNL

All cards are BaseCard objects.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from math import log, exp, ceil
import numpy as np
from numpy import unique, hstack

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.cards.base_card import BaseCard
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    string_or_blank, blank, fields, components_or_blank
)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16


class DELAY(BaseCard):
    type = 'DELAY'

    def __init__(self, sid, nodes, components, delays, comment=''):
        """
        +-------+-----+-----------+-----+--------+------+-----+--------+-----+
        |   1   |  2  |     3     |  4  |   5    |  6   |  7  |   8    |  9  |
        +=======+=====+===========+=====+========+======+=====+========+=====+
        | DELAY | SID | POINT ID1 | C1  |   T1   | P2   | C2  |   T2   |     |
        +-------+-----+-----------+-----+--------+------+-----+--------+-----+
        """
        if comment:
            self.comment = comment

        #: Identification number of DELAY entry. (Integer > 0)
        self.sid = sid
        #: Grid, extra, or scalar point identification number. (Integer > 0)
        self.nodes = nodes
        #: Component number. (Integers 1 through 6 for grid points; zero or blank for extra
        #: or scalar points)
        self.components = components
        #: Time delay (tau) for designated point Pi and component Ci. (Real)
        self.delays = delays

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DELAY card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nodes = [integer(card, 2, 'node')]
        components = [integer(card, 3, 'components')]
        delays = [double_or_blank(card, 4, 'delay')]
        assert components[0] in [0, 1, 2, 3, 4, 5, 6], components
        if card.field(5):
            nodes.append(integer(card, 5, 'node'))
            components.append(integer(card, 6, 'components'))
            delays.append(double_or_blank(card, 7, 'delay'))
            assert components[1] in [0, 1, 2, 3, 4, 5, 6], components
        return DELAY(sid, nodes, components, delays, comment=comment)

    def add(self, delay):
        assert self.sid == delay.sid, 'sid=%s delay.sid=%s' % (self.sid, delay.sid)
        if delay.comment:
            if hasattr('_comment'):
                self._comment += delay.comment
            else:
                self._comment = delay.comment
        self.nodes += delay.nodes
        self.components += delay.components
        self.delays += delay.delays

    def get_delay_at_freq(self, freq):
        return self.nodes, self.components, self.delays

    #def cross_reference(self, model):
        #"""
        #Cross links the card so referenced cards can be extracted directly

        #Parameters
        #----------
        #model : BDF()
            #the BDF object
        #"""
        #msg = ', which is required by DELAY sid=%s' % self.sid
        #self.nodes_ref = model.Node(self.node_ids, msg=msg)

    #@property
    #def node_id1(self):
        #if isinstance(self.nodes[0], integer_types):
            #return self.nodes[0]
        #return self.nodes_ref[0].nid

    #@property
    #def node_id2(self):
        #if isinstance(self.nodes[1], integer_types):
            #return self.nodes[1]
        #return self.nodes_ref[1].nid

    @property
    def node_ids(self):
        node_ids = [self.node_id1]
        if len(self.components) == 2:
            node_ids.append(self.node_id2)
        return node_ids

    def raw_fields(self):
        list_fields = ['DELAY', self.sid]
        for nid, comp, delay in zip(self.node_ids, self.components, self.delays):
            if isinstance(nid, integer_types):
                nidi = nid
            else:
                nidi = nid.nid
            list_fields += [nidi, comp, delay]
        return list_fields

    def write_card(self, size=8, is_double=False):
        msg = self.comment
        node_ids = self.node_ids
        if size == 8:
            for nid, comp, delay in zip(node_ids, self.components, self.delays):
                msg += print_card_8(['DELAY', self.sid, nid, comp, delay])
        else:
            for nid, comp, delay in zip(node_ids, self.components, self.delays):
                msg += print_card_16(['DELAY', self.sid, nid, comp, delay])
        return msg


class DPHASE(BaseCard):
    type = 'DPHASE'

    def __init__(self, sid, nodes, components, phase_leads, comment=''):
        """
        +--------+-----+-----------+-----+------+------+-----+-----+-----+
        |   1    |  2  |     3     |  4  |  5   |  6   |  7  |  8  |  9  |
        +========+=====+===========+=====+======+======+=====+=====+=====+
        | DPHASE | SID | POINT ID1 | C1  | TH1  |  P2  | C2  | TH2 |     |
        +--------+-----+-----------+-----+------+------+-----+-----+-----+
        """
        if comment:
            self.comment = comment
        self.sid = sid
        self.nodes = nodes
        self.components = components
        self.phase_leads = phase_leads

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a DPHASE card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nodes = [integer(card, 2, 'node')]
        components = [integer(card, 3, 'components')]
        phase_leads = [double_or_blank(card, 4, 'phase_lead')]
        assert components[0] in [0, 1, 2, 3, 4, 5, 6], components
        if card.field(5):
            nodes.append(integer(card, 5, 'node'))
            components.append(integer(card, 6, 'components'))
            phase_leads.append(double_or_blank(card, 7, 'phase_lead'))
            assert components[1] in [0, 1, 2, 3, 4, 5, 6], components
        return DPHASE(sid, nodes, components, phase_leads, comment=comment)

    def add(self, dphase):
        assert self.sid == dphase.sid, 'sid=%s dphase.sid=%s' % (self.sid, dphase.sid)
        if dphase.comment:
            if hasattr('_comment'):
                self._comment += dphase.comment
            else:
                self._comment = dphase.comment
        self.nodes += dphase.nodes
        self.components += dphase.components
        self.phase_leads += dphase.phase_leads

    #def cross_reference(self, model):
        #"""
        #Cross links the card so referenced cards can be extracted directly

        #Parameters
        #----------
        #model : BDF()
            #the BDF object
        #"""
        #msg = ', which is required by DPHASE sid=%s' % self.sid
        #self.nodes_ref = model.Nodes(self.node_ids, msg=msg)

    #@property
    #def node_id1(self):
        #if isinstance(self.nodes[0], integer_types):
            #return self.nodes[0]
        #return self.nodes_ref[0].nid

    #@property
    #def node_id2(self):
        #if isinstance(self.nodes[1], integer_types):
            #return self.nodes[1]
        #return self.nodes_ref[1].nid

    @property
    def node_ids(self):
        node_ids = [self.node_id1]
        if len(self.components) == 2:
            node_ids.append(self.node_id2)
        return node_ids

    def raw_fields(self):
        list_fields = ['DPHASE', self.sid]
        for nid, comp, delay in zip(self.nodes, self.components, self.phase_leads):
            if isinstance(nid, integer_types):
                nidi = nid
            else:
                nidi = nid.nid
            list_fields += [nidi, comp, delay]
        return list_fields

    def write_card(self, size=8, is_double=False):
        msg = self.comment
        node_ids = self.node_ids
        if size == 8:
            for nid, comp, delay in zip(node_ids, self.components, self.phase_leads):
                msg += print_card_8(['DPHASE', self.sid, nid, comp, delay])
        else:
            for nid, comp, delay in zip(node_ids, self.components, self.phase_leads):
                msg += print_card_16(['DPHASE', self.sid, nid, comp, delay])
        return msg


class FREQ(BaseCard):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems.

    +------+-----+-----+-----+------+-----+-----+-----+-----+
    |  1   |  2  |  3  |  4  |  5   |  6  |  7  |  8  |  9  |
    +======+=====+=====+=====+======+=====+=====+=====+=====+
    | FREQ | SID | F1  | F2  | etc. |     |     |     |     |
    +------+-----+-----+-----+------+-----+-----+-----+-----+
    """
    type = 'FREQ'

    def __init__(self, sid, freqs, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.freqs = np.unique(freqs)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FREQ card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        freqs = fields(double, card, 'freq', i=2, j=len(card))
        return FREQ(sid, freqs, comment=comment)

    def get_freqs(self):
        return self.freqs

    def add_frequencies(self, freqs):
        """
        Combines the frequencies from 1 FREQx object with another.
        All FREQi entries with the same frequency set identification numbers
        will be used. Duplicate frequencies will be ignored.

        Parameters
        ----------
        freqs : ???
            the frequencies for a FREQx object
        """
        #print("self.freqs = ",self.freqs)
        #print("freqs = ",freqs)
        self.freqs = unique(hstack([self.freqs, freqs]))

    def add_frequency_object(self, freq):
        """
        :param freq: a FREQx object

        .. seealso:: :func:`addFrequencies`
        """
        self.add_frequencies(freq.freqs)

    def raw_fields(self):
        list_fields = ['FREQ', self.sid] + list(self.freqs)
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class FREQ1(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, frequency
    increment, and the number of increments desired.

    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |   1   |  2  | 3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+=====+
    | FREQ1 | SID |  F1 | DF  | NDF |     |     |     |     |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    """
    type = 'FREQ1'

    def __init__(self, sid, f1, df, ndf, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.f1 = f1
        self.df = df
        self.ndf = ndf

        freqs = []
        for i in range(ndf):
            freqs.append(f1 + i * df)
        self.freqs = unique(freqs)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FREQ1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        f1 = double_or_blank(card, 2, 'f1', 0.0)
        df = double(card, 3, 'df')
        ndf = integer_or_blank(card, 4, 'ndf', 1)
        assert len(card) <= 5, 'len(FREQ card) = %i\ncard=%s' % (len(card), card)
        return FREQ1(sid, f1, df, ndf, comment=comment)

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class FREQ2(FREQ):
    """
    Defines a set of frequencies to be used in the solution of frequency
    response problems by specification of a starting frequency, final
    frequency, and the number of logarithmic increments desired.

    +-------+-----+-----+-----+-----+-----+-----+-----+-----+
    |   1   |  2  | 3   |  4  |  5  |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+=====+=====+=====+=====+=====+
    | FREQ2 | SID |  F1 | F2  | NDF |     |     |     |     |
    +-------+-----+-----+-----+-----+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    """
    type = 'FREQ2'

    def __init__(self, sid, f1, f2, ndf=1, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.f1 = f1
        self.f2 = f2
        self.ndf = ndf

        d = 1. / ndf * log(f2 / f1)
        freqs = []
        for i in range(ndf):
            freqs.append(f1 * exp(i * d))  # 0 based index
        self.freqs = np.unique(freqs)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FREQ2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        f1 = double(card, 2, 'f1')  # default=0.0 ?
        f2 = double(card, 3, 'f2')
        ndf = integer_or_blank(card, 4, 'nf', 1)
        assert len(card) <= 5, 'len(FREQ2 card) = %i\ncard=%s' % (len(card), card)
        return FREQ2(sid, f1, f2, ndf, comment=comment)
        #return FREQ(sid, freqs, comment=comment)


class FREQ3(FREQ):
    """
    +-------+-----+------+-------+--------+-----+---------+
    |   1   |  2  |   3  |   4   |    5   |  6  |    7    |
    +=======+=====+======+=======+========+=====+=========+
    | FREQ3 | SID |  F1  |  F2   |  TYPE  | NEF | CLUSTER |
    +-------+-----+------+-------+--------+-----+---------+
    | FREQ3 |  6  | 20.0 | 200.0 | LINEAR | 10  |   2.0   |
    +-------+-----+------+-------+--------+-----+---------+
    """
    type = 'FREQ3'

    def __init__(self, f1, f2=None, Type='LINEAR', nef=10, cluster=1.0, comment=''):
        if comment:
            self.comment = comment
        if f2 is None:
            f2 = f1
        self.sid = sid
        self.f1 = f1
        self.f2 = f2
        self.Type = Type
        self.nef = nef
        self.cluster = cluster

    @classmethod
    def add_card(cls, card, comment=''):
        sid = integer(card, 1, 'sid')
        f1 = double(card, 1, 'f1')
        f2 = integer_or_blank(card, 1, 'f2', f1)
        Type = string_or_blank(card, 1, 'Type', 'LINEAR')
        nef = integer_or_blank(card, 1, 'nef', 10)
        cluster = double_or_blank(card, 1, 'cluster', 1.0)

        return FREQ3(sid, f1, f2, Type, nef, cluster, comment='')

    def raw_fields(self):
        return ['FREQ3', self.sid, self.f1, self.f2, self.Type, self.nef, self.cluster]

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class FREQ4(FREQ):
    """
    Defines a set of frequencies used in the solution of modal frequency
    response problems by specifying the amount of 'spread' around each natural
    frequency and the number of equally spaced excitation frequencies within
    the spread.

    +-------+-----+-----+-----+------+-----+-----+-----+-----+
    |   1   |  2  | 3   |  4  |  5   |  6  |  7  |  8  |  9  |
    +=======+=====+=====+=====+======+=====+=====+=====+=====+
    | FREQ4 | SID |  F1 | F2  | FSPD | NFM |     |     |     |
    +-------+-----+-----+-----+------+-----+-----+-----+-----+

    .. note:: this card rewrites as a FREQ card
    .. todo:: not done...
    """
    type = 'FREQ4'

    def __init__(self, sid, f1, f2, fspread, nfm, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.f1 = f1
        self.f2 = f2
        self.fspread = fspread
        self.nfm = nfm

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a FREQ4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        f1 = double_or_blank(card, 2, 'f1', 0.0)
        f2 = double_or_blank(card, 3, 'f2', 1.e20)
        fspread = double_or_blank(card, 4, 'fspd', 0.1)
        nfm = integer_or_blank(card, 5, 'nfm', 3)
        assert len(card) <= 6, 'len(FREQ card) = %i\ncard=%s' % (len(card), card)
        return FREQ4(sid, f1, f2, fspread, nfm, comment=comment)

    def raw_fields(self):
        list_fields = ['FREQ4', self.sid, self.f1, self.f2, self.fspread,
                       self.nfm]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


#class FREQ5(FREQ):
    #type = 'FREQ5'

    #def __init__(self, card=None, data=None, comment=''):
        #if comment:
            # self.comment = comment
        #raise NotImplementedError()

    #def write_card(self, size=8, is_double=False):
        #card = self.repr_fields()
        #if size == 8:
            #return self.comment + print_card_8(card)
        #return self.comment + print_card_16(card)


class NLPARM(BaseCard):
    """
    Defines a set of parameters for nonlinear static analysis iteration
    strategy.

    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |    1   |  2     |  3   |  4   |    5    |   6   |   7     |  8      |   9    |
    +========+========+======+======+=========+=======+=========+=========+========+
    | NLPARM |   ID   | NINC |  DT  | KMETHOD | KSTEP | MAXITER |  CONV   | INTOUT |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |        |  ESPU  | EPSP | EPSW | MAXDIV  | MAXQN | MAXLS   | FSTRESS | LSTOL  |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    |        | MAXBIS |      |      |         | MAXR  |         | RTOLB   | CONV   |
    +--------+--------+------+------+---------+-------+---------+---------+--------+
    """
    type = 'NLPARM'

    def __init__(self, nlparm_id, ninc=10, dt=0.0, kmethod='AUTO', kstep=5,
                 max_iter=25, conv='PW', int_out='NO',
                 eps_u=0.01, eps_p=0.01, eps_w=0.01, max_div=3, max_qn=None, max_ls=4,
                 fstress=0.2, ls_tol=0.5, max_bisect=5, max_r=20., rtol_b=20., comment=''):
        if comment:
            self.comment = comment
        self.nlparm_id = nlparm_id
        self.ninc = ninc
        self.dt = dt
        self.kmethod = kmethod
        self.kstep = kstep
        self.max_iter = max_iter
        self.conv = conv
        self.int_out = int_out

        # line 2
        self.eps_p = eps_p
        self.eps_u = eps_u
        self.eps_w = eps_w
        self.max_div = max_div
        self.max_qn = max_qn
        self.max_ls = max_ls
        self.fstress = fstress
        self.ls_tol = ls_tol

        # line 3
        self.max_bisect = max_bisect
        self.max_r = max_r
        self.rtol_b = rtol_b

        if self.max_qn is None:
            if kmethod == 'PFNT':
                self.max_qn = 0
            else:
                self.max_qn = max_iter

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NLPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nlparm_id = integer(card, 1, 'nlparm_id')
        ninc = integer_or_blank(card, 2, 'ninc', 10)
        dt = double_or_blank(card, 3, 'dt', 0.0)
        kmethod = string_or_blank(card, 4, 'kmethod', 'AUTO')
        kstep = integer_or_blank(card, 5, 'kstep', 5)
        max_iter = integer_or_blank(card, 6, 'max_iter', 25)
        conv = string_or_blank(card, 7, 'conv', 'PW')
        int_out = string_or_blank(card, 8, 'intOut', 'NO')

        # line 2
        eps_u = double_or_blank(card, 9, 'eps_u', 0.01)
        eps_p = double_or_blank(card, 10, 'eps_p', 0.01)
        eps_w = double_or_blank(card, 11, 'eps_w', 0.01)
        max_div = integer_or_blank(card, 12, 'max_div', 3)

        if kmethod == 'PFNT':
            max_qn = integer_or_blank(card, 13, 'max_qn', 0)
        else:
            max_qn = integer_or_blank(card, 13, 'max_qn', max_iter)

        max_ls = integer_or_blank(card, 14, 'max_ls', 4)
        fstress = double_or_blank(card, 15, 'fstress', 0.2)
        ls_tol = double_or_blank(card, 16, 'ls_tol', 0.5)

        # line 3
        max_bisect = integer_or_blank(card, 17, 'max_bisect', 5)
        max_r = double_or_blank(card, 21, 'max_r', 20.)
        rtol_b = double_or_blank(card, 23, 'rtol_b', 20.)
        assert len(card) <= 24, 'len(NLPARM card) = %i\ncard=%s' % (len(card), card)
        return NLPARM(nlparm_id, ninc, dt, kmethod, kstep, max_iter, conv,
                      int_out, eps_u, eps_p, eps_w, max_div,
                      max_qn, max_ls, fstress,
                      ls_tol, max_bisect, max_r,
                      rtol_b, comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a NLPARM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (nlparm_id, ninc, dt, kmethod, kstep, max_iter, conv, int_out, eps_u, eps_p,
         eps_w, max_div, max_qn, max_ls, fstress, ls_tol, max_bisect, max_r,
         rtol_b) = data

        if kmethod == 1:
            kmethod = 'AUTO'
        elif kmethod == 2:
            kmethod = 'ITER'
        elif kmethod == 4:
            kmethod = 'SEMI'
        elif kmethod == 3:
            kmethod = 'ADAPT'
        else:
            msg = 'nlparm_id=%s kmethod=%r data=%s' % (nlparm_id, kmethod, data)
            raise NotImplementedError(msg)

        if conv == 1:
            conv = 'W'
        elif conv == 2:
            conv = 'P'
        elif conv == 3:
            conv = 'PW'
        elif conv == 4:
            conv = 'U'
        elif conv == 5:
            conv = 'UW'
        elif conv == 6:
            conv = 'UP'
        elif conv == 7:
            conv = 'UPW'
        else:
            msg = 'nlparm_id=%s conv=%r data=%s' % (nlparm_id, conv, data)
            raise NotImplementedError(msg)

        if int_out == 0:
            int_out = 'NO'
        elif int_out == 1:
            int_out = 'YES'
        elif int_out == 2:
            int_out = 'ALL'
        else:
            msg = 'nlparm_id=%s int_out=%r data=%s' % (nlparm_id, int_out, data)
            raise NotImplementedError(msg)
        return NLPARM(nlparm_id, ninc, dt, kmethod, kstep, max_iter, conv,
                      int_out, eps_u, eps_p, eps_w, max_div,
                      max_qn, max_ls, fstress,
                      ls_tol, max_bisect, max_r,
                      rtol_b, comment=comment)

    def raw_fields(self):
        list_fields = ['NLPARM', self.nlparm_id, self.ninc, self.dt, self.kmethod,
                       self.kstep, self.max_iter, self.conv, self.int_out, self.eps_u,
                       self.eps_p, self.eps_w, self.max_div, self.max_qn, self.max_ls,
                       self.fstress, self.ls_tol, self.max_bisect, None, None, None,
                       self.max_r, None, self.rtol_b]
        return list_fields

    def repr_fields(self):
        ninc = set_blank_if_default(self.ninc, 10)
        dt = set_blank_if_default(self.dt, 0.0)
        kmethod = set_blank_if_default(self.kmethod, 'AUTO')
        kstep = set_blank_if_default(self.kstep, 5)
        max_iter = set_blank_if_default(self.max_iter, 25)
        conv = set_blank_if_default(self.conv, 'PW')
        int_out = set_blank_if_default(self.int_out, 'NO')
        eps_u = set_blank_if_default(self.eps_u, 0.01)
        eps_p = set_blank_if_default(self.eps_p, 0.01)
        eps_w = set_blank_if_default(self.eps_w, 0.01)
        max_div = set_blank_if_default(self.max_div, 3)
        max_qn = set_blank_if_default(self.max_qn, self.max_iter)
        max_ls = set_blank_if_default(self.max_ls, 4)
        fstress = set_blank_if_default(self.fstress, 0.2)
        ls_tol = set_blank_if_default(self.ls_tol, 0.5)
        max_bisect = set_blank_if_default(self.max_bisect, 5)
        max_r = set_blank_if_default(self.max_r, 20.)
        rtol_b = set_blank_if_default(self.rtol_b, 20.)

        list_fields = ['NLPARM', self.nlparm_id, ninc, dt, kmethod, kstep, max_iter,
                       conv, int_out, eps_u, eps_p, eps_w, max_div, max_qn, max_ls,
                       fstress, ls_tol, max_bisect, None, None, None, max_r, None,
                       rtol_b]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card) # having trouble with double precision...
        return self.comment + print_card_16(card)


class NLPCI(BaseCard):
    type = 'NLPCI'

    def __init__(self, nlpci_id, Type='CRIS', minalr=0.25, maxalr=4.,
                 scale=0., desiter=12, mxinc=20, comment=''):
        if comment:
            self.comment = comment
        self.nlpci_id = nlpci_id
        self.Type = Type
        self.minalr = minalr
        self.maxalr = maxalr
        self.scale = scale
        self.desiter = desiter
        self.mxinc = mxinc

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a NLPCI card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nlpci_id = integer(card, 1, 'nlpci_id')
        Type = string_or_blank(card, 2, 'Type', 'CRIS')
        minalr = double_or_blank(card, 3, 'minalr', 0.25)
        maxalr = double_or_blank(card, 4, 'maxalr', 4.0)
        scale = double_or_blank(card, 5, 'scale', 0.0)
        blank(card, 6, 'blank')
        desiter = integer_or_blank(card, 7, 'desiter', 12)
        mxinc = integer_or_blank(card, 8, 'mxinc', 20)
        return NLPCI(nlpci_id, Type=Type, minalr=minalr, maxalr=maxalr,
                     scale=scale, desiter=desiter, mxinc=mxinc, comment=comment)

    def raw_fields(self):
        list_fields = ['NLPCI', self.nlpci_id, self.Type, self.minalr,
                       self.maxalr, self.scale, None, self.desiter, self.mxinc]
        return list_fields

    def repr_fields(self):
        #minalr = set_blank_if_default(self.minalr, 0.25)
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class TF(BaseCard):
    """
    Defines a dynamic transfer function of the form:
        (B0 + B1 p + B2 *p2)*ud  sum(A0_i + A1_i*p + A2_i*p2)*ui = 0

    +----+-----+-----+------+------+------+--------+----+----+
    |  1 |  2  |  3  |   4  |   5  |   6  |    7   |  8 |  9 |
    +====+=====+=====+======+======+======+========+====+====+
    | TF | SID | GD  |  CD  |  B0  |  B1  |   B2   |    |    |
    +----+-----+-----+------+------+------+--------+----+----+
    |    | G_1 | C_1 | A0_1 | A1_1 | A2_1 | -etc.- |    |    |
    +----+-----+-----+------+------+------+--------+----+----+

    """
    type = 'TF'
    def __init__(self, sid, nid0, c, b0, b1, b2, nids, components, a, comment=''):
        if comment:
            self.comment = comment
        self.sid = sid
        self.nid0 = nid0
        self.c = c
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.nids = nids
        self.components = components
        self.a = a

    def validate(self):
        pass
        #assert len(self.grids1) > 0, 'ngrids1=%s\n%s' % (len(self.grids1), str(self))

    #def cross_reference(self, model):
        #pass

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        nid0 = integer(card, 2, 'nid0')
        # component 0 means an SPOINT/EPOINT
        c = components_or_blank(card, 3, 'components_0', 0)
        b0 = double_or_blank(card, 4, 'b0', 0.)
        b1 = double_or_blank(card, 5, 'b1', 0.)
        b2 = double_or_blank(card, 6, 'b2', 0.)

        nfields = len(card) - 9
        nrows = nfields // 8
        if nfields % 8 > 0:
            nrows += 1

        nids = []
        components = []
        a = []
        for irow in range(nrows):
            j = irow * 8 + 9
            #ifield = irow + 1
            nid = integer(card, j, 'grid_%i' % (irow + 1))
            component = components_or_blank(card, j + 1, 'components_%i' % (irow + 1), 0)
            a0 = double_or_blank(card, j + 2, 'a0_%i' % (irow + 1), 0.)
            a1 = double_or_blank(card, j + 3, 'a1_%i' % (irow + 1), 0.)
            a2 = double_or_blank(card, j + 4, 'a2_%i' % (irow + 1), 0.)
            nids.append(nid)
            components.append(component)
            a.append([a0, a1, a2])
        return TF(sid, nid0, c, b0, b1, b2, nids, components, a,
                  comment=comment)

    def raw_fields(self):
        list_fields = ['TF', self.sid, self.nid0, self.c, self.b0, self.b1, self.b2, None, None]
        for grid, c, (a0, a1, a2) in zip(self.nids, self.components, self.a):
            list_fields += [grid, c, a0, a1, a2, None, None, None]
        return list_fields

    def write_card(self, size=8, is_double=False):
        # double precision?
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class TSTEP(BaseCard):
    """
    Transient Time Step
    Defines time step intervals at which a solution will be generated and
    output in transient analysis.

    +-------+------+------+------+------+-----+-----+-----+-----+
    |   1   |   2  |  3   |  4   |  5   |  6  |  7  |  8  |  9  |
    +=======+======+======+======+======+=====+=====+=====+=====+
    | TSTEP | SID  |  N1  | DT1  | NO1  |     |     |     |     |
    +-------+------+------+------+------+-----+-----+-----+-----+
    |       |      |  N2  | DT2  | NO2  |     |     |     |     |
    +-------+------+------+------+------+-----+-----+-----+-----+
    |       |      | etc. |      |      |     |     |     |     |
    +-------+------+------+------+------+-----+-----+-----+-----+

    +-------+------+------+------+------+-----+-----+-----+-----+
    |   1   |   2  |  3   |  4   |  5   |  6  |  7  |  8  |  9  |
    +=======+======+======+======+======+=====+=====+=====+=====+
    | TSTEP | 101  | 9000 | .001 | 9000 |     |     |     |     |
    +-------+------+------+------+------+-----+-----+-----+-----+
    |       |      | 1000 | .001 | 1    |     |     |     |     |
    +-------+------+------+------+------+-----+-----+-----+-----+
    """
    type = 'TSTEP'

    def __init__(self, sid, N, DT, NO, comment=''):
        """
        Creates a TSTEP card

        Parameters
        ----------
        sid : int
            the time step id
        N : List[int/None]
            ???
        DT : List[float/None]
            ???
        NO : List[int/None]
            ???
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment
        self.sid = sid
        #: Number of time steps of value DTi. (Integer > 1)
        self.N = N
        #: Time increment (float)
        self.DT = DT
        #: Skip factor for output. Every NOi-th step will be saved for output (default=1)
        self.NO = NO

    def validate(self):
        assert len(self.N) == len(self.DT), 'N=%s DT=%s' % (self.N, self.DT)
        assert len(self.N) == len(self.NO), 'N=%s NO=%s' % (self.N, self.NO)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TSTEP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        N = []
        DT = []
        NO = []

        nrows = int(ceil((len(card) - 1.) / 8.))
        for i in range(nrows):
            n = 8 * i + 1
            ni = integer_or_blank(card, n + 1, 'N' + str(i), 1)
            dt = double_or_blank(card, n + 2, 'dt' + str(i), 0.)
            no = integer_or_blank(card, n + 3, 'NO' + str(i), 1)
            N.append(ni)
            DT.append(dt)
            NO.append(no)
        return TSTEP(sid, N, DT, NO, comment=comment)

    def raw_fields(self):
        list_fields = ['TSTEP', self.sid]
        for (N, dt, no) in zip(self.N, self.DT, self.NO):
            list_fields += [N, dt, no, None, None, None, None, None]
        return list_fields

    def repr_fields(self):
        return self.raw_fields()

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)


class TSTEPNL(BaseCard):
    """
    Defines parametric controls and data for nonlinear transient structural or
    heat transfer analysis. TSTEPNL is intended for SOLs 129, 159, and 600.
    Parameters for Nonlinear Transient Analysis.

    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |    1    |   2    |   3    |   4   |   5    |   6    |   7   |    8    |  9   |
    +=========+========+========+=======+========+========+=======+=========+======+
    | TSTEPNL |   ID   |  NDT   |  DT   |   NO   | METHOD | KSTEP | MAXITER | CONV |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |         |  ESPU  |  EPSP  |  EPSW | MAXDIV | MAXQN  | MAXLS | FSTRESS |      |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+
    |         | MAXBIS | ADJUST | MSTEP |   RB   |  MAXR  | UTOL  | RTOLB   |      |
    +---------+--------+--------+-------+--------+--------+-------+---------+------+

    method = None for NX, but apparently TSTEP as well, which is not in the QRG
    """
    type = 'TSTEPNL'
    allowed_methods = ['AUTO', 'ITER', 'ADAPT', 'SEMI', 'FNT', 'PFNT', # MSC
                       'TSTEP'] # NX

    def __init__(self, sid, ndt, dt, no, method='ADAPT', kstep=None,
                 max_iter=10, conv='PW', eps_u=1.e-2, eps_p=1.e-3,
                 eps_w=1.e-6, max_div=2, max_qn=10, max_ls=2,
                 fstress=0.2, max_bisect=5, adjust=5, mstep=None,
                 rb=0.6, max_r=32., utol=0.1, rtol_b=20.,
                 min_iter=None, comment=''):
        """
        Creates a TSTEPNL card

        Parameters
        ----------
        sid : int
            the time step id
        ndt : ???
            ???
        dt : ???
            ???
        no : ???
            ???
        eps_u : float; default=1.e-2
            ???
        eps_p : float; default=1.e-3
            ???
        eps_w : float; default=1.e-6
            ???
        max_div : int; default=2
            ???
        max_qn : int; default=10
            ???
        max_ls : int; default=2
            ???
        fstress : float; default=0.2
            ???
        max_bisect : int; default=5
            ???
        adjust : int; default=5
            ???
        mstep : int; default=None
            ???
        rb : float; default=0.6
            ???
        max_r = float; default=32.
            ???
        utol = float; default=0.1
            ???
        rtol_b = float; default=20.
            ???
        min_iter : int; default=None
            not listed in all QRGs
        comment : str; default=''
            a comment for the card
        """
        if comment:
            self.comment = comment

        # line 1
        self.sid = sid
        self.ndt = ndt
        self.dt = dt
        self.no = no
        self.method = method
        self.kstep = kstep
        self.max_iter = max_iter
        self.conv = conv

        self.eps_u = eps_u
        self.eps_p = eps_p
        self.eps_w = eps_w
        self.max_div = max_div
        self.max_qn = max_qn
        self.max_ls = max_ls
        self.fstress = fstress

        # line 3
        self.max_bisect = max_bisect
        self.adjust = adjust
        self.mstep = mstep
        self.rb = rb
        self.max_r = max_r
        self.utol = utol
        self.rtol_b = rtol_b
        self.min_iter = min_iter
        assert self.ndt >= 3
        assert self.dt > 0.

    def validate(self):
        if self.method not in self.allowed_methods:
            msg = 'method=%r allowed_methods=[%s]' % (
                self.method, ', '.join(self.allowed_methods))
            raise ValueError(msg)

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a TSTEPNL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        sid = integer(card, 1, 'sid')
        ndt = integer(card, 2, 'ndt')
        dt = double(card, 3, 'dt')
        no = integer_or_blank(card, 4, 'no', 1)

        #: .. note:: not listed in all QRGs
        method = string_or_blank(card, 5, 'method', 'ADAPT')
        if method == 'ADAPT':
            kstep = integer_or_blank(card, 6, 'kStep', 2)
        elif method == 'ITER':
            kstep = integer_or_blank(card, 6, 'kStep', 10)
        elif method in ['AUTO', 'TSTEP', 'SEMI']:
            kstep = None
            #kstep = blank(card, 6, 'kStep') #: .. todo:: not blank
        else:
            msg = 'invalid TSTEPNL Method.  method=%r; allowed_methods=[%s]' % (
                method, ', '.join(cls.allowed_methods))
            raise RuntimeError(msg)
        max_iter = integer_or_blank(card, 7, 'maxIter', 10)
        conv = string_or_blank(card, 8, 'conv', 'PW')

        # line 2
        eps_u = double_or_blank(card, 9, 'epsU', 1.E-2)
        eps_p = double_or_blank(card, 10, 'epsP', 1.E-3)
        eps_w = double_or_blank(card, 11, 'epsW', 1.E-6)
        max_div = integer_or_blank(card, 12, 'maxDiv', 2)
        max_qn = integer_or_blank(card, 13, 'maxQn', 10)
        max_ls = integer_or_blank(card, 14, 'MaxLs', 2)
        fstress = double_or_blank(card, 15, 'fStress', 0.2)

        # line 3
        max_bisect = integer_or_blank(card, 17, 'maxBisect', 5)
        adjust = integer_or_blank(card, 18, 'adjust', 5)
        mstep = integer_or_blank(card, 19, 'mStep')
        rb = double_or_blank(card, 20, 'rb', 0.6)
        max_r = double_or_blank(card, 21, 'maxR', 32.)
        utol = double_or_blank(card, 22, 'uTol', 0.1)
        rtol_b = double_or_blank(card, 23, 'rTolB', 20.)

        # not listed in all QRGs
        min_iter = integer_or_blank(card, 24, 'minIter')
        assert len(card) <= 25, 'len(TSTEPNL card) = %i\ncard=%s' % (len(card), card)
        return TSTEPNL(
            sid, ndt, dt, no, method, kstep, max_iter, conv,
            eps_u, eps_p, eps_w, max_div, max_qn, max_ls, fstress,
            max_bisect, adjust, mstep, rb, max_r, utol, rtol_b, min_iter,
            comment=comment)

    @classmethod
    def add_op2_data(cls, data, comment=''):
        """
        Adds a TSTEPNL card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        comment : str; default=''
            a comment for the card
        """
        (sid, ndt, dt, no, method, kstep, max_iter, conv, eps_u, eps_p, eps_w,
         max_div, max_qn, max_ls, fstress, max_bisect,
         adjust, mstep, rb, max_r, utol, rtol_b) = data

        if method == 1:
            method = 'AUTO'
        elif method == 3:
            method = 'ADAPT'
        else:
            raise NotImplementedError('tstepnl=%s method=%r data=%s' % (sid, method, data))

        if conv == 3:
            conv = 'PW'
        elif conv == 4:
            conv = 'U'
        #elif conv == 3:
            #conv = 'ADAPT'
        else:
            raise NotImplementedError('tstepnl=%s conv=%r data=%s' % (sid, conv, data))

        min_iter = None  # not listed in DMAP 2005
        return TSTEPNL(
            sid, ndt, dt, no, method, kstep, max_iter, conv,
            eps_u, eps_p, eps_w, max_div, max_qn, max_ls, fstress,
            max_bisect, adjust, mstep, rb, max_r, utol, rtol_b, min_iter,
            comment=comment)
        #self.sid = sid
        #self.ndt = ndt
        #self.dt = dt
        #self.no = no
        #self.method = method
        #self.kStep = kStep
        #self.maxIter = maxIter
        #self.conv = conv

        ## line 2
        #self.epsU = epsU
        #self.epsP = epsP
        #self.epsW = epsW
        #self.maxDiv = maxDiv
        #self.maxQn = maxQn
        #self.MaxLs = maxLs
        #self.fStress = fStress

        ## line 3
        #self.maxBisect = maxBisect
        #self.adjust = adjust
        #self.mStep = mStep
        #self.rb = rb
        #self.maxR = maxR
        #self.uTol = uTol
        #self.rTolB = rTolB

    def raw_fields(self):
        list_fields = ['TSTEPNL', self.sid, self.ndt, self.dt, self.no,
                       self.method, self.kstep, self.max_iter, self.conv, self.eps_u,
                       self.eps_p, self.eps_w, self.max_div, self.max_qn, self.max_ls,
                       self.fstress, None, self.max_bisect, self.adjust, self.mstep,
                       self.rb, self.max_r, self.utol, self.rtol_b, self.min_iter]
        return list_fields

    def repr_fields(self):
        #no = set_blank_if_default(self.no,1)
        no = self.no
        method = set_blank_if_default(self.method, 'ADAPT')

        kstep = self.kstep
        #if self.method == 'ADAPT':
            #kStep = set_blank_if_default(self.kStep, 2)
        #elif self.method == 'ITER':
            #kStep = set_blank_if_default(self.kStep, 10)
        #else:
            #msg = 'invalid TSTEPNL Method.  method=|%s|' %(self.method)
            #raise RuntimeError(msg)

        #maxIter = set_blank_if_default(self.maxIter, 10)
        conv = set_blank_if_default(self.conv, 'PW')

        eps_u = set_blank_if_default(self.eps_u, 1e-2)
        eps_p = set_blank_if_default(self.eps_p, 1e-3)
        eps_w = set_blank_if_default(self.eps_w, 1e-6)
        max_div = set_blank_if_default(self.max_div, 2)
        max_qn = set_blank_if_default(self.max_qn, 10)
        max_ls = set_blank_if_default(self.max_ls, 2)
        fstress = set_blank_if_default(self.fstress, 0.2)

        max_bisect = set_blank_if_default(self.max_bisect, 5)
        adjust = set_blank_if_default(self.adjust, 5)
        rb = set_blank_if_default(self.rb, 0.6)
        max_r = set_blank_if_default(self.max_r, 32.)
        utol = set_blank_if_default(self.utol, 0.1)
        rtol_b = set_blank_if_default(self.rtol_b, 20.)

        list_fields = ['TSTEPNL', self.sid, self.ndt, self.dt, no, method,
                       kstep, self.max_iter, conv, eps_u, eps_p, eps_w, max_div, max_qn,
                       max_ls, fstress, None, max_bisect, adjust, self.mstep, rb,
                       max_r, utol, rtol_b, self.min_iter]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
