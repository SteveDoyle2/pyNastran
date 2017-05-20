# pylint: disable=R0904,R0902,E1101,E1103,C0111,C0302,C0103,W0101
from __future__ import print_function
from six import string_types
import numpy as np
from numpy.linalg import norm

from pyNastran.utils import integer_types
from pyNastran.bdf.cards.elements.bars import CBAR, LineElement
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_string_or_blank)
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

class CBEAM(CBAR):
    """

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | X1  | X2  | X3  | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    or

    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |   1   |  2  |  3  |  4  |  5  |  6  |  7  |  8  |    9     |
    +=======+=====+=====+=====+=====+=====+=====+=====+==========+
    | CBEAM | EID | PID | GA  | GB  | G0  |     |     | OFFT/BIT |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | PA  | PB  | W1A | W2A | W3A | W1B | W2B | W3B      |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+
    |       | SA  | SB  |     |     |     |     |     |          |
    +-------+-----+-----+-----+-----+-----+-----+-----+----------+

    offt/bit are MSC specific fields
    """
    type = 'CBEAM'
    _field_map = {
        1: 'eid', 2:'pid', 3:'ga', 4:'gb', #5:'x_g0', 6:'g1', 7:'g2',
        #8:'offt',
        9:'pa', 10:'pb',
        17:'sa', 18:'sb',
    }

    def _update_field_helper(self, n, value):
        if n == 11:
            self.wa[0] = value
        elif n == 12:
            self.wa[1] = value
        elif n == 13:
            self.wa[2] = value

        elif n == 14:
            self.wb[0] = value
        elif n == 15:
            self.wb[1] = value
        elif n == 16:
            self.wb[2] = value
        else:
            if self.g0 is not None:
                if n == 5:
                    self.g0 = value
                else:  # offt
                    msg = 'Field %r=%r is an invalid %s entry or is unsupported.' % (
                        n, value, self.type)
                    raise KeyError(msg)
            else:
                if n == 5:
                    self.x[0] = value
                elif n == 6:
                    self.x[1] = value
                elif n == 7:
                    self.x[2] = value
                else:
                    msg = 'Field %r=%r is an invalid %s entry or is unsupported.' % (
                        n, value, self.type)
                    raise KeyError(msg)

    def __init__(self, eid, pid, nids, x, g0, offt='GGG', bit=None,
                 pa=0, pb=0, wa=None, wb=None, sa=0, sb=0, comment=''):
        """
        Adds a CBEAM card

        Parameters
        ----------
        pid : int
            property id
        mid : int
            material id
        nids : List[int, int]
            node ids; connected grid points at ends A and B
        x : List[float, float, float]
            Components of orientation vector, from GA, in the displacement
            coordinate system at GA (default), or in the basic coordinate system
        g0 : int
            Alternate method to supply the orientation vector using grid
            point G0. Direction of is from GA to G0. is then transferred
            to End A
        offt : str; default='GGG'
            Offset vector interpretation flag
            None : bit is active
        bit : float; default=None
            Built-in twist of the cross-sectional axes about the beam axis
            at end B relative to end A.
            For beam p-elements ONLY!
            None : offt is active
        pa / pb : int; default=0
            Pin Flag at End A/B.  Releases the specified DOFs
        wa / wb : List[float, float, float]
            Components of offset vectors from the grid points to the end
            points of the axis of the shear center
        sa / sb : int; default=0
            Scalar or grid point identification numbers for the ends A and B,
            respectively. The degrees-of-freedom at these points are the
            warping variables . SA and SB cannot be specified for
            beam p-elements
        comment : str; default=''
            a comment for the card

        offt/bit are MSC specific fields
        """
        LineElement.__init__(self)
        if comment:
            self.comment = comment
        if wa is None:
            wa = np.zeros(3, dtype='float64')
        else:
            wa = np.asarray(wa)
        if wb is None:
            wb = np.zeros(3, dtype='float64')
        else:
            wb = np.asarray(wb)

        self.eid = eid
        self.pid = pid
        self.ga = nids[0]
        self.gb = nids[1]
        self.x = x
        self.g0 = g0
        self.offt = offt
        self.bit = bit
        self.pa = pa
        self.pb = pb
        self.wa = wa
        self.wb = wb
        self.sa = sa
        self.sb = sb
        self._validate_input()

    @classmethod
    def add_card(cls, card, comment=''):
        """
        Adds a CBEAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        eid = integer(card, 1, 'eid')
        pid = integer_or_blank(card, 2, 'pid', eid)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')

        x, g0 = cls._init_x_g0(card, eid)
        offt, bit = cls._init_offt_bit(card, eid)# offt doesn't exist in NX nastran
        pa = integer_or_blank(card, 9, 'pa', 0)
        pb = integer_or_blank(card, 10, 'pb', 0)

        wa = np.array([double_or_blank(card, 11, 'w1a', 0.0),
                       double_or_blank(card, 12, 'w2a', 0.0),
                       double_or_blank(card, 13, 'w3a', 0.0)], 'float64')

        wb = np.array([double_or_blank(card, 14, 'w1b', 0.0),
                       double_or_blank(card, 15, 'w2b', 0.0),
                       double_or_blank(card, 16, 'w3b', 0.0)], 'float64')

        sa = integer_or_blank(card, 17, 'sa', 0)
        sb = integer_or_blank(card, 18, 'sb', 0)
        assert len(card) <= 19, 'len(CBEAM card) = %i\ncard=%s' % (len(card), card)
        return CBEAM(eid, pid, [ga, gb], x, g0, offt, bit,
                     pa=pa, pb=pb, wa=wa, wb=wb, sa=sa, sb=sb, comment=comment)

    @classmethod
    def add_op2_data(cls, data, f, comment=''):
        """
        Adds a CBEAM card from the OP2

        Parameters
        ----------
        data : List[varies]
            a list of fields defined in OP2 format
        f : int
            beam flag
            0 : basic
                [x1, x2, x3] is used
            1 : cid
                [x1, x2, x3] is used
            2 : grid
                g0 is used instead of [x1, x2, x3]
        comment : str; default=''
            a comment for the card
        """
        #: .. todo:: verify
        assert len(data) == 2, 'data=%s len(data)=%s' % (data, len(data))
        #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
        #        [f,g0]]
        #data = [[eid,pid,ga,gb,sa,sb, pa,pb,w1a,w2a,w3a,w1b,w2b,w3b],
        #        [f,x1,x2,x3]]

        main, aft = data
        flag = aft[0]
        assert f == flag, 'f=%s flag=%s' % (f, flag)
        if flag == 0:
            # basic cid
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
            assert len(aft) == 4, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            x1, x2, x3 = aft[1:]
            g0 = None
            x = np.array([x1, x2, x3], dtype='float64')
        elif flag == 1:
            # global cid
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, x1, x2, x3]]
            assert len(aft) == 4, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            g0 = None
            x1, x2, x3 = aft[1:]
            x = np.array([x1, x2, x3], dtype='float64')
        elif flag == 2:
            # grid option
            #data_in = [[eid, pid, ga, gb, sa, sb, pa, pb, w1a, w2a, w3a, w1b, w2b, w3b],
                       #[f, g0]]
            assert len(aft) == 2, 'f=%s aft=%s len(aft)=%s' % (f, aft, len(aft))
            g0 = data[1][1]
            x = None
        else:
            raise NotImplementedError()

        eid = main[0]
        pid = main[1]
        ga = main[2]
        gb = main[3]
        sa = main[4]
        sb = main[5]

        #offt = str(data[6]) # GGG
        bit = None # ???
        offt = 'GGG'  #: .. todo:: is this correct???

        pa = main[6]
        pb = main[7]

        wa = np.array([main[8], main[9], main[10]], 'float64')
        wb = np.array([main[11], main[12], main[13]], 'float64')
        return CBEAM(eid, pid, [ga, gb], x, g0, offt, bit,
                     pa=pa, pb=pb, wa=wa, wb=wb, sa=sa, sb=sb, comment=comment)

    def _validate_input(self):
        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

    def Nodes(self):
        return [self.ga, self.gb]

    @classmethod
    def _init_offt_bit(cls, card, eid):
        """
        offt doesn't exist in NX nastran
        """
        field8 = integer_double_string_or_blank(card, 8, 'field8')
        if isinstance(field8, float):
            offt = None
            bit = field8
        elif field8 is None:
            offt = 'GGG'  # default
            bit = None
        elif isinstance(field8, string_types):
            bit = None
            offt = field8
            msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
            assert offt[0] in ['G', 'B', 'O', 'E'], msg
            assert offt[1] in ['G', 'B', 'O', 'E'], msg
            assert offt[2] in ['G', 'B', 'O', 'E'], msg
        else:
            msg = ('field8 on %s card is not a string(offt) or bit '
                   '(float)...field8=%s\n' % (cls.type, field8))
            raise RuntimeError("Card Instantiation: %s" % msg)
        return offt, bit

    def Centroid(self):
        return (self.ga_ref.get_position() + self.gb_ref.get_position()) / 2.

    def center_of_mass(self):
        """
        A          B
        *----------*
        ^          ^
        | wa       | wb
        |          |
        1----------2

        1-2 are the nodes of the bar
        A-B defines the axis of the shear center


             ^ z
        +--+ |
        |  | |
        |  | 2---> y
        |  +--+
        |     |
        +-----+
        """
        cda = self.ga_ref.cid_ref
        cdb = self.gb_ref.cid_ref
        ga = self.ga_ref.get_position() + cda.transform_node_to_global_assuming_rectangular(self.wa)
        gb = self.gb_ref.get_position() + cdb.transform_node_to_global_assuming_rectangular(self.wb)
        #x = self.get_orientation_vector()
        return (ga + gb) / 2.

    def get_axes(self, model, debug=False):
        """
        OFFT flag
        ---------
        ABC or A-B-C (an example is G-G-G or B-G-G)
        while the slots are:
         - A -> orientation; values=[G, B]
         - B -> End A; values=[G, O]
         - C -> End B; values=[G, O]

        and the values for A,B,C mean:
         - B -> basic
         - G -> global
         - O -> orientation

        so for example G-G-G, that's global for all terms.
        BOG means basic orientation, orientation end A, global end B

        so now we're left with what does basic/global/orientation mean?
        - basic -> the global coordinate system defined by cid=0
        - global -> the local coordinate system defined by the
                    CD field on the GRID card, but referenced by
                    the CBAR/CBEAM
        - orientation -> ???

        NX Nastran uses GGG implicitly
        """
        eid = self.eid
        (nid1, nid2) = self.node_ids
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]
        n1 = node1.get_position()
        n2 = node2.get_position()
        centroid = (n1 + n2) / 2.
        i = n2 - n1
        Li = norm(i)
        ihat = i / Li

        is_failed = True
        if self.g0:
            #debug = False
            msg = 'which is required by %s eid=%s\n%s' % (self.type, self.g0, str(self))
            g0_ref = model.Node(self.g0, msg=msg)
            if debug:  # pragma: no cover
                print('  g0 = %s' % self.g0)
                print('  g0_ref = %s' % g0_ref)
            n0 = g0_ref.get_position()
            v = n0 - n1
        else:
            #debug = False
            ga = model.nodes[self.Ga()]
            cda = ga.Cd()
            cda_ref = model.Coord(cda)
            v = cda_ref.transform_node_to_global(self.x)
            if debug:  # pragma: no cover
                print('  ga = %s' % self.ga)
                if cda != 0:
                    print('  cd = %s' % cda_ref)
                else:
                    print('  cd = 0')

                print('  x = %s' % self.x)
                print('  v = %s' % v)
            #v = self.x

        offt_vector, offt_end_a, offt_end_b = self.offt
        if debug:  # pragma: no cover
            print('  offt vector,A,B=%r' % (self.offt))
        # if offt_end_a == 'G' or (offt_end_a == 'O' and offt_vector == 'G'):

        cd1 = node1.Cd()
        cd2 = node2.Cd()
        cd1_ref = model.Coord(cd1)
        cd2_ref = model.Coord(cd2)
        # node1.cd_ref, node2.cd_ref

        if offt_vector == 'G':
            # end A
            # global - cid != 0
            if cd1 != 0:
                v = cd1_ref.transform_node_to_global_assuming_rectangular(v)
                #if node1.cd_ref.type not in ['CORD2R', 'CORD1R']:
                    #msg = 'invalid Cd type (%r) on Node %i; expected CORDxR' % (
                        #node1.cd_ref.type, node1.nid)
                    #self.log.error(msg)
                    #continue
                    #raise NotImplementedError(node1.cd)
        elif offt_vector == 'B':
            # basic - cid = 0
            pass
        else:
            msg = 'offt_vector=%r is not supported; offt=%s' % (offt_vector, self.offt)
            #self.log.error(msg)
            return is_failed, msg
            #raise NotImplementedError(msg)
        #print('v = %s' % v)

        # rotate wa
        wa = self.wa
        if offt_end_a == 'G':
            if cd1 != 0:
                #if node1.cd.type not in ['CORD2R', 'CORD1R']:
                    #continue # TODO: support CD transform
                # TODO: fixme
                wa = cd1_ref.transform_node_to_global_assuming_rectangular(wa)
        elif offt_end_a == 'B':
            pass
        elif offt_end_a == 'O':
            # TODO: fixme
            wa = cd1_ref.transform_node_to_global_assuming_rectangular(n1 - wa)
        else:
            msg = 'offt_end_a=%r is not supported; offt=%s' % (offt_end_a, self.offt)
            self.log.error(msg)
            return is_failed, msg
            #raise NotImplementedError(msg)

        #print('wa = %s' % wa)
        # rotate wb
        wb = self.wb
        if offt_end_b == 'G':
            if cd2 != 0:
                #if cd2_ref.type not in ['CORD2R', 'CORD1R']:
                    #continue # TODO: MasterModelTaxi
                # TODO: fixme
                wb = cd2_ref.transform_node_to_global_assuming_rectangular(wb)

        elif offt_end_b == 'B':
            pass
        elif offt_end_b == 'O':
            # TODO: fixme
            wb = cd1_ref.transform_node_to_global_assuming_rectangular(n2 - wb)
        else:
            msg = 'offt_end_b=%r is not supported; offt=%s' % (offt_end_b, self.offt)
            model.log.error(msg)
            return is_failed, msg
            #raise NotImplementedError(msg)

        #print('wb =', wb)
        ## concept has a GOO
        #if not self.offt in ['GGG', 'BGG']:
            #msg = 'offt=%r for CBAR/CBEAM eid=%s is not supported...skipping' % (
                #self.offt, eid)
            #self.log.error(msg)
            #continue

        vhat = v / norm(v) # j
        try:
            z = np.cross(ihat, vhat) # k
        except ValueError:
            msg = 'Invalid vector length\n'
            msg += 'n1  =%s\n' % str(n1)
            msg += 'n2  =%s\n' % str(n2)
            msg += 'nid1=%s\n' % str(nid1)
            msg += 'nid2=%s\n' % str(nid2)
            msg += 'i   =%s\n' % str(i)
            msg += 'Li  =%s\n' % str(Li)
            msg += 'ihat=%s\n' % str(ihat)
            msg += 'v   =%s\n' % str(v)
            msg += 'vhat=%s\n' % str(vhat)
            msg += 'z=cross(ihat, vhat)'
            print(msg)
            raise ValueError(msg)

        zhat = z / norm(z)
        yhat = np.cross(zhat, ihat) # j
        if debug:
            print('  centroid = %s' % centroid)
            print('  ihat = %s' % ihat)
            print('  yhat = %s' % yhat)
            print('  zhat = %s' % zhat)
        #if eid == 5570:
            #print('  check - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                  #eid, yhat, zhat, v, i, nid1, n1, nid2, n2))

        if norm(ihat) == 0.0 or norm(yhat) == 0.0 or norm(z) == 0.0:
            print('  invalid_orientation - eid=%s yhat=%s zhat=%s v=%s i=%s n%s=%s n%s=%s' % (
                eid, yhat, zhat, v, i, nid1, n1, nid2, n2))
        elif not np.allclose(norm(yhat), 1.0) or not np.allclose(norm(zhat), 1.0) or Li == 0.0:
            print('  length_error        - eid=%s Li=%s Lyhat=%s Lzhat=%s'
                  ' v=%s i=%s n%s=%s n%s=%s' % (
                      eid, Li, norm(yhat), norm(zhat), v, i, nid1, n1, nid2, n2))

        #print('adding bar %s' % bar_type)
        #print('   centroid=%s' % centroid)
        #print('   yhat=%s len=%s' % (yhat, np.linalg.norm(yhat)))
        #print('   zhat=%s len=%s' % (zhat, np.linalg.norm(zhat)))
        #print('   Li=%s' % (Li))
        is_failed = False
        return is_failed, wa, wb, ihat, yhat, zhat

    def Mid(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Mid()

    def Area(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Area()

    def Nsm(self):
        if isinstance(self.pid, integer_types):
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Nsm()

    @property
    def is_offt(self):
        """is the offt flag active?"""
        if isinstance(self.offt, string_types):
            return True
        assert isinstance(self.bit, float), 'bit=%s type=%s' % (self.bit, type(self.bit))
        return False

    @property
    def is_bit(self):
        """is the bit flag active?"""
        return not self.is_offt

    def get_offt_bit_defaults(self):
        """
        offt doesn't exist in NX nastran
        """
        if self.is_offt:
            field8 = set_blank_if_default(self.offt, 'GGG')
        else:
            field8 = set_blank_if_default(self.bit, 0.0)
        return field8

    def cross_reference(self, model):
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object
        """
        msg = ' which is required by CBEAM eid=%s' % (self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.ga_ref = self.ga
        self.gb = model.Node(self.gb, msg=msg)
        self.gb_ref = self.gb
        self.nodes = model.Nodes([self.ga.nid, self.gb.nid], msg=msg)
        self.nodes_ref = self.nodes
        self.pid = model.Property(self.pid, msg=msg)
        self.pid_ref = self.pid
        if self.g0:
            g0 = model.nodes[self.g0]
            self.g0_vector = g0.get_position() - self.ga.get_position()
        else:
            self.g0_vector = self.x
        if model.is_nx:
            assert self.offt == 'GGG', 'NX only support offt=GGG; offt=%r' % self.offt

    def safe_cross_reference(self, model):
        msg = ' which is required by CBEAM eid=%s' % (self.eid)
        self.ga = model.Node(self.ga, msg=msg)
        self.gb = model.Node(self.gb, msg=msg)

        self.ga_ref = self.ga
        self.gb_ref = self.gb
        try:
            self.pid = model.Property(self.pid, msg=msg)
            self.pid_ref = self.pid
        except KeyError:
            model.log.warning('pid=%s%s' % (self.pid, msg))

        if self.g0:
            try:
                g0 = model.nodes[self.g0]
                self.g0_vector = g0.get_position() - self.ga.get_position()
            except KeyError:
                model.log.warning('Node=%s%s' % (self.g0, msg))
        else:
            self.g0_vector = self.x

    def uncross_reference(self):
        self.pid = self.Pid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        del self.ga_ref, self.gb_ref, self.pid_ref

    def _verify(self, xref=False):
        eid = self.eid
        pid = self.Pid()
        edges = self.get_edge_ids()
        if xref:  # True
            mid = self.Mid()
            nsm = self.Nsm()
            assert isinstance(mid, int), 'mid=%r' % mid
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert self.pid_ref.type in ['PBEAM', 'PBEAML', 'PBCOMP'], '%s%s' % (self, self.pid_ref)
            A = self.Area()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
            assert isinstance(A, float), 'eid=%s A=%r' % (eid, A)
            assert isinstance(L, float), 'eid=%s L=%r' % (eid, L)
            assert isinstance(mpl, float), 'eid=%s mass_per_length=%r' % (eid, mpl)
            assert isinstance(mass, float), 'eid=%s mass=%r' % (eid, mass)
            assert L > 0.0, 'eid=%s L=%s' % (eid, L)

    def raw_fields(self):
        (x1, x2, x3) = self.get_x_g0_defaults()
        offt = self.get_offt_bit_defaults()
        ga, gb = self.node_ids
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                       self.pa, self.pb] + list(self.wa) + list(self.wb) + [self.sa, self.sb]
        return list_fields

    def repr_fields(self):
        w1a = set_blank_if_default(self.wa[0], 0.0)
        w2a = set_blank_if_default(self.wa[1], 0.0)
        w3a = set_blank_if_default(self.wa[2], 0.0)
        w1b = set_blank_if_default(self.wb[0], 0.0)
        w2b = set_blank_if_default(self.wb[1], 0.0)
        w3b = set_blank_if_default(self.wb[2], 0.0)
        pa = set_blank_if_default(self.pa, 0)
        pb = set_blank_if_default(self.pb, 0)

        sa = set_blank_if_default(self.sa, 0)
        sb = set_blank_if_default(self.sb, 0)
        (x1, x2, x3) = self.get_x_g0_defaults()
        offt = self.get_offt_bit_defaults()
        ga, gb = self.node_ids
        list_fields = ['CBEAM', self.eid, self.Pid(), ga, gb, x1, x2, x3, offt,
                       pa, pb, w1a, w2a, w3a,
                       w1b, w2b, w3b, sa, sb]
        return list_fields

    def write_card(self, size=8, is_double=False):
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)
