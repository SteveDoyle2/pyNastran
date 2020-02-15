# pylint: disable=C0103
"""
defines:
 - CBEAM
 - BEAMOR

"""
from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from numpy.linalg import norm  # type: ignore

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.elements.bars import (
    LineElement, init_x_g0, BaseCard, rotate_v_wa_wb, check_offt)
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank, integer_double_string_or_blank,
    integer_double_or_blank, integer_string_or_blank,
)
from pyNastran.bdf.field_writer_8 import set_blank_if_default
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.utils.mathematics import integrate_positive_unit_line
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF


class CBEAM(LineElement):
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

    @classmethod
    def export_to_hdf5(cls, h5_file, model, eids):
        """exports the elements in a vectorized way"""
        encoding = model._encoding
        #comments = []
        pids = []
        nodes = []

        x = []
        g0 = []
        offt = []
        bit = []

        pa = []
        pb = []
        wa = []
        wb = []
        sa = []
        sb = []
        nan = np.full(3, np.nan)
        for eid in eids:
            element = model.elements[eid]
            #comments.append(element.comment)
            pids.append(element.pid)
            nodes.append(element.nodes)
            if element.g0 is None:
                x.append(element.x)
                g0.append(-1)
            else:
                x.append(nan)
                g0.append(element.g0)

            if element.bit is not None:
                bit.append(element.bit)
                offt.append(b'')
            else:
                bit.append(np.nan)
                offti = element.offt
                if isinstance(offti, integer_types):
                    offti = str(offti)
                offt.append(offti.encode(encoding))


            pa.append(element.pa)
            pb.append(element.pb)
            sa.append(element.sa)
            sb.append(element.sb)
            wa.append(element.wa)
            wb.append(element.wb)
        #h5_file.create_dataset('_comment', data=comments)
        h5_file.create_dataset('eid', data=eids)
        h5_file.create_dataset('nodes', data=nodes)
        h5_file.create_dataset('pid', data=pids)
        #print('x =', x)
        #print('g0 =', g0)
        h5_file.create_dataset('x', data=x)
        h5_file.create_dataset('g0', data=g0)
        h5_file.create_dataset('offt', data=offt)
        h5_file.create_dataset('bit', data=bit)

        h5_file.create_dataset('pa', data=pa)
        h5_file.create_dataset('pb', data=pb)

        h5_file.create_dataset('sa', data=sa)
        h5_file.create_dataset('sb', data=sb)

        h5_file.create_dataset('wa', data=wa)
        h5_file.create_dataset('wb', data=wb)


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

        if isinstance(offt, str):
            offt = offt.replace('E', 'O')
            offt = int(offt) if offt.isdigit() else offt
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
        self.ga_ref = None
        self.gb_ref = None
        self.pid_ref = None
        self.g0_ref = None
        self.g0_vector = None

    def validate(self):
        msg = ''
        if self.x is None:
            if not isinstance(self.g0, integer_types):
                msg += 'CBEAM eid=%s: x is None, so g0=%s must be an integer' % (self.eid, self.g0)
        else:
            if not isinstance(self.x, (list, np.ndarray)):
                msg += 'CBEAM eid=%s: x=%s and g0=%s, so x must be a list; type(x)=%s' % (
                    self.eid, self.x, self.g0, type(self.x))
        if msg:
            raise ValueError(msg)

        if self.g0 is not None:
            assert isinstance(self.g0, integer_types), 'g0=%s must be an integer' % self.g0
        if self.g0 in [self.ga, self.gb]:
            msg = 'G0=%s cannot be GA=%s or GB=%s' % (self.g0, self.ga, self.gb)
            raise RuntimeError(msg)

        if self.bit is None and self.offt is None:
            msg = 'OFFT/BIT must not be None; offt=%r bit=%s' % (self.offt, self.bit)
            raise RuntimeError(msg)

        if self.offt is not None:
            if isinstance(self.offt, integer_types):
                assert self.offt in [1, 2, 21, 22, 41, 42], 'invalid offt; offt=%i' % self.offt
                #raise NotImplementedError('invalid offt; offt=%i' % self.offt)
            elif isinstance(self.offt, str):
                check_offt(self)
            else:
                raise TypeError('invalid offt expected a string of length 3 '
                                'offt=%r; Type=%s' % (self.offt, type(self.offt)))

    @classmethod
    def add_card(cls, card, beamor=None, comment=''):
        """
        Adds a CBEAM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        beamor : BEAMOR() or None
            defines the defaults
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')

        pid_default = eid
        x1_default, x2_default, x3_default = 0., 0., 0.
        offt_default = 'GGG'
        if beamor is not None:
            if beamor.pid is not None:
                pid_default = beamor.pid
            if beamor.x is None:
                x1_default = beamor.g0
                x2_default = None
                x3_default = None
            else:
                x1_default, x2_default, x3_default = beamor.x
            offt_default = beamor.offt

        pid = integer_or_blank(card, 2, 'pid', pid_default)
        ga = integer(card, 3, 'ga')
        gb = integer(card, 4, 'gb')

        x, g0 = init_x_g0(card, eid, x1_default, x2_default, x3_default)
        offt, bit = _init_offt_bit(card, eid, offt_default)# offt doesn't exist in NX nastran
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

    def Nodes(self):
        return [self.ga, self.gb]

    def Centroid(self):
        """"""
        node1 = self.ga_ref
        node2 = self.gb_ref
        xyz1 = node1.get_position()
        xyz2 = node2.get_position()
        centroid = (xyz1 + xyz2) / 2.
        return centroid

    def center_of_mass(self):
        """the centroid formuala is way more complicated if you consider the nonstructural mass axis"""
        elem = self
        prop = self.pid_ref
        node1 = self.ga_ref
        node2 = self.gb_ref
        xyz1 = node1.get_position()
        xyz2 = node2.get_position()
        #centroid = ( + self.gb_ref.get_position()) / 2.
        centroid = (xyz1 + xyz2) / 2.
        #length = norm(xyz2 - xyz1)
        #cda = model.nodes[n1].cid_ref
        #cdb = model.nodes[n2].cid_ref

        model = None
        log = None
        is_failed, out = elem.get_axes_by_nodes(model, self.pid_ref, node1, node2, xyz1, xyz2, log)
        if is_failed:
            #model.log.error(out)
            raise RuntimeError(out)

        wa, wb, _ihat, jhat, khat = out
        p1 = xyz1 + wa
        p2 = xyz2 + wb

        if prop.type == 'PBEAM':
            rho = prop.Rho()

            # we don't call the MassPerLength method so we can put the NSM centroid
            # on a different axis (the PBEAM is weird)
            mass_per_lengths = []
            nsm_per_lengths = []
            for (area, nsm) in zip(prop.A, prop.nsm):
                mass_per_lengths.append(area * rho)
                nsm_per_lengths.append(nsm)
            mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
            nsm_per_length = integrate_positive_unit_line(prop.xxb, nsm_per_lengths)
            nsm_n1 = (p1 + jhat * prop.m1a + khat * prop.m2a)
            nsm_n2 = (p2 + jhat * prop.m1b + khat * prop.m2b)
            #print("nsm_per_length=%s" % nsm_per_length)
            #print("nsm_n1=%s" % nsm_n1)
            #print("nsm_n2=%s" % nsm_n2)
            nsm_centroid = (nsm_n1 + nsm_n2) / 2.
            #if nsm != 0.:
                #p1_nsm = p1 + prop.ma
                #p2_nsm = p2 + prop.mb
        elif prop.type == 'PBEAML':
            mass_per_lengths = prop.get_mass_per_lengths()
            #mass_per_length = prop.MassPerLength() # includes simplified nsm

            # m1a, m1b, m2a, m2b=0.
            nsm_centroid = (p1 + p2) / 2.

            # mass_per_length already includes nsm
            mass_per_length = integrate_positive_unit_line(prop.xxb, mass_per_lengths)
            nsm_per_length = 0.

            #print('mass_per_lengths=%s nsm_per_lengths=%s' % (
                #mass_per_lengths, nsm_per_lengths))
            #print('mass_per_length=%s nsm_per_length=%s' % (
                #mass_per_length, nsm_per_length))

            #nsm_centroid = np.zeros(3) # TODO: what is this...
            #nsm = prop.nsm[0] * length # TODO: simplified
        elif prop.type == 'PBCOMP':
            mass_per_length = prop.MassPerLength()
            nsm_per_length = prop.nsm
            nsm_n1 = (p1 + jhat * prop.m1 + khat * prop.m2)
            nsm_n2 = (p2 + jhat * prop.m1 + khat * prop.m2)
            nsm_centroid = (nsm_n1 + nsm_n2) / 2.
        #elif prop.type == 'PBMSECT':
            #continue
            #mass_per_length = prop.MassPerLength()
            #m = mass_per_length * length
            #nsm = prop.nsm
        elif prop.type == 'PBMSECT':
            mass_per_length = 0.  ## TODO: fix me
            nsm_per_length = prop.nsm
            nsm_centroid = (p1 + p2) / 2.
        else:
            raise NotImplementedError(prop.type)

        total_mass = mass_per_length + nsm_per_length
        if total_mass == 0.0:
            return centroid
        centroid2 = (centroid * mass_per_length + nsm_centroid * nsm_per_length) / total_mass
        return centroid2

    def center_of_mass_xform(self):
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

    def get_axes(self, model):
        """
        Gets the axes of a CBAR/CBEAM, while respecting the OFFT flag.

        Notes
        -----
        :func:`pyNastran.bdf.cards.elements.bars.rotate_v_wa_wb` for a
        description of the OFFT flag.

        is_passed: bool
        out: (wa, wb, ihat, jhat, khat)
        """
        is_failed = True
        #TODO: not integrated with CBAR yet...
        if self.bit is not None:
            print(self.get_stats())
            return is_failed, None

        check_offt(self)
        is_failed = True
        ihat = None
        yhat = None
        zhat = None

        eid = self.eid
        (nid1, nid2) = self.node_ids
        node1 = model.nodes[nid1]
        node2 = model.nodes[nid2]
        xyz1 = node1.get_position()
        xyz2 = node2.get_position()

        elem = model.elements[eid]
        pid_ref = elem.pid_ref
        if pid_ref is None:
            pid_ref = model.Property(elem.pid)
        assert not isinstance(pid_ref, integer_types), elem

        is_failed, (wa, wb, ihat, yhat, zhat) = self.get_axes_by_nodes(
            model, pid_ref, node1, node2, xyz1, xyz2, model.log)
        return is_failed, (wa, wb, ihat, yhat, zhat)

    def get_axes_by_nodes(self, model, pid_ref, node1, node2, xyz1, xyz2, log):
        """
        Gets the axes of a CBAR/CBEAM, while respecting the OFFT flag.

        Notes
        -----
        :func:`pyNastran.bdf.cards.elements.bars.rotate_v_wa_wb` for a
        description of the OFFT flag.

        """
        #TODO: not integrated with CBAR yet...

        is_failed = True
        eid = self.eid
        #centroid = (n1 + n2) / 2.
        #i = n2 - n1
        #Li = norm(i)
        #ihat = i / Li

        elem = self
        #(nid1, nid2) = elem.node_ids
        #node1 = model.nodes[nid1]
        #node2 = model.nodes[nid2]
        #xyz1 = node1.get_position()
        #xyz2 = node2.get_position()

        # wa/wb are not considered in i_offset
        # they are considered in ihat
        i = xyz2 - xyz1
        Li = norm(i)
        if Li == 0.:
            msg = 'xyz1=%s xyz2=%s\n%s' % (xyz1, xyz2, self)
            raise ValueError(msg)
        i_offset = i / Li

        unused_v, wa, wb, xform = rotate_v_wa_wb(
            model, elem,
            xyz1, xyz2, node1, node2,
            i_offset, i, eid, Li, log)
        if wb is None:
            # one or more of v, wa, wb are bad

            # xform is xform_offset...assuming None
            ihat = None
            yhat = None
            zhat = None
            return is_failed, (wa, wb, ihat, yhat, zhat)

        ihat = xform[0, :]
        yhat = xform[1, :]
        zhat = xform[2, :]

        is_failed = False
        return is_failed, (wa, wb, ihat, yhat, zhat)

    @property
    def node_ids(self):
        return [self.Ga(), self.Gb()]

    def get_edge_ids(self):
        return [tuple(sorted(self.node_ids))]

    @property
    def nodes(self):
        return [self.ga, self.gb]

    @nodes.setter
    def nodes(self, values):
        self.ga = values[0]
        self.gb = values[1]

    @property
    def nodes_ref(self):
        return [self.ga_ref, self.gb_ref]

    @nodes_ref.setter
    def nodes_ref(self, values):
        assert values is not None, values
        self.ga_ref = values[0]
        self.gb_ref = values[1]

    def Mid(self):
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Mid()

    def Area(self):
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Area()

    def Nsm(self):
        if self.pid_ref is None:
            raise RuntimeError('Element eid=%i has not been '
                               'cross referenced.\n%s' % (self.eid, str(self)))
        return self.pid_ref.Nsm()

    @property
    def is_offt(self):
        """is the offt flag active?"""
        if self.bit is not None:
            assert isinstance(self.bit, float), 'bit=%r type=%s' % (self.bit, type(self.bit))
            return False
        #assert isinstance(self.offt, str), 'offt=%r' % self.offt
        return True

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

    def cross_reference(self, model: BDF) -> None:
        """
        Cross links the card so referenced cards can be extracted directly

        Parameters
        ----------
        model : BDF()
            the BDF object

        """
        msg = ', which is required by CBEAM eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.nodes_ref = [self.ga_ref, self.gb_ref]
        self.pid_ref = model.Property(self.pid, msg=msg)
        if self.g0:
            self.g0_ref = model.nodes[self.g0]
            self.g0_vector = self.g0_ref.get_position() - self.ga_ref.get_position()
        else:
            self.g0_vector = self.x
        if model.is_nx:
            assert self.offt == 'GGG', 'NX only support offt=GGG; offt=%r' % self.offt

    def safe_cross_reference(self, model, xref_errors):
        msg = ', which is required by CBEAM eid=%s' % (self.eid)
        self.ga_ref = model.Node(self.ga, msg=msg)
        self.gb_ref = model.Node(self.gb, msg=msg)
        self.nodes_ref = [self.ga_ref, self.gb_ref]
        self.pid_ref = model.safe_property(self.pid, self.eid, xref_errors, msg=msg)

        if self.g0:
            try:
                self.g0_ref = model.nodes[self.g0]
                self.g0_vector = self.g0_ref.get_position() - self.ga_ref.get_position()
            except KeyError:
                model.log.warning('Node=%s%s' % (self.g0, msg))
        else:
            self.g0_vector = self.x

    def uncross_reference(self) -> None:
        """Removes cross-reference links"""
        self.pid = self.Pid()
        self.ga = self.Ga()
        self.gb = self.Gb()
        self.g0 = self.G0()
        self.ga_ref = None
        self.gb_ref = None
        self.g0_ref = None
        self.pid_ref = None

    def _verify(self, xref):
        eid = self.eid
        unused_pid = self.Pid()
        unused_edges = self.get_edge_ids()
        if xref:  # True
            prop = self.pid_ref
            assert prop.type in ['PBEAM', 'PBEAML', 'PBCOMP', 'PBMSECT'], prop
            mid = self.Mid()
            nsm = self.Nsm()
            assert isinstance(mid, int), 'mid=%r' % mid
            assert isinstance(nsm, float), 'nsm=%r' % nsm
            assert self.pid_ref.type in ['PBEAM', 'PBEAML', 'PBCOMP', 'PBMSECT'], '%s%s' % (self, self.pid_ref)
            A = self.Area()
            mpl = self.MassPerLength()
            L = self.Length()
            mass = self.Mass()
            assert isinstance(A, float), 'eid=%s A=%r' % (eid, A)
            assert isinstance(L, float), 'eid=%s L=%r' % (eid, L)
            assert isinstance(mpl, float), 'eid=%s mass_per_length=%r' % (eid, mpl)
            assert isinstance(mass, float), 'eid=%s mass=%r' % (eid, mass)
            assert L > 0.0, 'eid=%s L=%s' % (eid, L)

    def Ga(self):
        """gets Ga/G1"""
        if self.ga_ref is None:
            return self.ga
        return self.ga_ref.nid

    def Gb(self):
        """gets Gb/G2"""
        if self.gb_ref is None:
            return self.gb
        return self.gb_ref.nid

    def G0(self):
        """gets G0"""
        if self.g0_ref is None:
            return self.g0
        return self.g0_ref.nid

    def get_x_g0_defaults(self):
        """
        X and G0 compete for the same fields, so the method exists to
        make it easier to write the card

        Returns
        -------
        x_g0 : varies
            g0 : List[int, None, None]
            x : List[float, float, float]

        Notes
        -----
        Used by CBAR and CBEAM

        """
        if self.g0 is not None:
            return (self.G0(), None, None)
        else:
            return list(self.x)

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

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)

    def write_card_16(self, is_double=False):
        card = self.repr_fields()
        return self.comment + print_card_16(card)

def _init_offt_bit(card, unused_eid, offt_default):
    """
    offt doesn't exist in NX nastran
    """
    field8 = integer_double_string_or_blank(card, 8, 'field8', offt_default)
    if isinstance(field8, float):
        offt = None
        bit = field8
    elif field8 is None:
        offt = 'GGG'  # default
        bit = None
    elif isinstance(field8, integer_types):
        bit = None
        offt = field8
    elif isinstance(field8, str):
        bit = None
        offt = field8
        msg = 'invalid offt parameter of CBEAM...offt=%s' % offt
        assert offt[0] in ['G', 'B', 'O', 'E'], msg
        assert offt[1] in ['G', 'B', 'O', 'E'], msg
        assert offt[2] in ['G', 'B', 'O', 'E'], msg
    else:
        msg = ('field8 on %s card is not a string(offt) or bit '
               '(float)...field8=%s\n' % (card.field(0), field8))
        raise SyntaxError("Card Instantiation: %s" % msg)
    return offt, bit


class BEAMOR(BaseCard):
    """
    +--------+-----+---+---+---+-------+-----+-------+------+
    |    1   |  2  | 3 | 4 | 5 |   6   |  7  |   8   |  9   |
    +========+=====+===+===+===+=======+=====+=======+======+
    | BEAMOR | PID |   |   |   | G0/X1 |  X2 |  X3   | OFFT |
    +--------+-----+---+---+---+-------+-----+-------+------+
    | BEAMOR | 39  |   |   |   |  0.6  | 2.9 | -5.87 | GOG  |
    +--------+-----+---+---+---+-------+-----+-------+------+

    """
    type = 'BEAMOR'
    def __init__(self, pid, is_g0, g0, x, offt='GGG', comment=''):
        BaseCard.__init__(self)
        if comment:
            self.comment = comment
        self.pid = pid
        self.g0 = g0
        self.x = x
        self.offt = offt

    @classmethod
    def _init_from_empty(cls):
        pid = 1
        is_g0 = True
        g0 = 1
        x = None
        return BEAMOR(pid, is_g0, g0, x, offt='GGG', comment='')

    @classmethod
    def add_card(cls, card, comment=''):
        pid = integer_or_blank(card, 2, 'pid')

        # x / g0
        field5 = integer_double_or_blank(card, 5, 'g0_x1', 0.0)
        if isinstance(field5, integer_types):
            is_g0 = True
            g0 = field5
            x = [0., 0., 0.]
        elif isinstance(field5, float):
            is_g0 = False
            g0 = None
            x = np.array([field5,
                          double_or_blank(card, 6, 'x2', 0.0),
                          double_or_blank(card, 7, 'x3', 0.0)],
                         dtype='float64')
        else:
            raise NotImplementedError('BEAMOR field5 = %r' % field5)
        offt = integer_string_or_blank(card, 8, 'offt', 'GGG')
        assert len(card) <= 9, 'len(BEAMOR card) = %i\ncard=%s' % (len(card), card)
        return BEAMOR(pid, is_g0, g0, x, offt=offt, comment=comment)

    def raw_fields(self):
        return ['BEAMOR', None, self.pid, None, None] + list(self.x) + [self.offt]

    def write_card(self, size: int=8, is_double: bool=False) -> str:
        card = self.repr_fields()
        if size == 8:
            return self.comment + print_card_8(card)
        return self.comment + print_card_16(card)
