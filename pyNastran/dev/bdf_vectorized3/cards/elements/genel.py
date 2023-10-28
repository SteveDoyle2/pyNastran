from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double,
    #integer_or_blank, double_or_blank, string_or_blank,
    #integer_double_or_blank,
    fields)
#from pyNastran.bdf.field_writer_8 import print_card_8 # , print_float_8, print_field_8

from pyNastran.dev.bdf_vectorized3.cards.base_card import Element, parse_element_check
from pyNastran.dev.bdf_vectorized3.cards.base_card import hslice_by_idim, make_idim # , VectorizedBaseCard, searchsorted_filter
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, get_print_card_size # , array_default_int


if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class GENEL(Element):
    """
    +-------+------+-----+------+------+------+------+-------+------+
    |   1   |   2  |  3  |   4  |   5  |   6  |   7  |   8   |   9  |
    +=======+======+=====+======+======+======+======+=======+======+
    | GENEL | EID  |     | UI1  | CI1  | UI2  | CI2  |  UI3  | CI3  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | UI4  | CI4 | UI5  | CI5  | etc. |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | UD   |     | UD1  | CD1  | UD2  | CD2  | etc.  |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | K/Z  |     | KZ11 | KZ21 | KZ31 | etc. | KZ22  | KZ32 |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | etc. |     | KZ33 | KZ43 | etc. |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       | S    |     | S11  | S12  | etc. |  S21 |  etc. |      |
    +-------+------+-----+------+------+------+------+-------+------+

    +-------+------+-----+------+------+------+------+-------+------+
    | GENEL |  629 |     |  1   |  1   |  13  |  4   |   42  |   0  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  24  |  2  |      |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  UD  |     |  6   |  2   |  33  |  0   |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  Z   | 1.0 | 2.0  | 3.0  | 4.0  | 5.0  |  6.0  | 7.0  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  8.0 | 9.0 | 10.0 |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  S   | 1.5 | 2.5  | 3.5  | 4.5  | 5.5  |  6.5  | 7.5  |
    +-------+------+-----+------+------+------+------+-------+------+
    |       |  8.5 |     |      |      |      |      |       |      |
    +-------+------+-----+------+------+------+------+-------+------+

    """
    type = 'GENEL'
    #pid = 0
    #_properties = ['node_ids', 'ul_nodes', 'ud_nodes', 'nodes']
    #@classmethod
    #def _init_from_empty(cls):
        #eid = 1
        #ul = None
        #ud = None
        #k = [1.]
        #z = None
        #return GENEL(eid, ul, ud, k, z, s=None, comment='')

    def __init__(self, model: BDF):
        super().__init__(model)
        #self.property_id = np.array([], dtype='int32')

    #def __init__(self, eid, ul, ud, k, z, s=None, comment=''):
        #"""creates a GENEL card

        #The required input is the {UL} list and the lower triangular
        #portion of [K] or [Z].  Additional input may include the {UD}
        #list and [S].  If [S] is input, must also be input.  If {UD} is
        #input but [S] is omitted, [S] is internally calculated. In this
        #case, {UD} must contain six and only six degrees-of freedom.
        #"""
        #BaseCard.__init__(self)
        #if comment:
            #self.comment = comment
        #self.eid = eid
        #self.ul = ul
        #self.ud = ud
        #if k is not None:
            #self.k = np.asarray(k)
            #self.z = None
        #else:
            #self.z = np.asarray(z)
            #self.k = None

        #if s is not None:
            #s = np.asarray(s)
        #self.s = s

        #self.ul_nodes_ref = None
        #self.ud_nodes_ref = None

    #def _finalize_hdf5(self, encoding):
        #self.ul = np.array(self.ul, dtype='int32')#.reshape(len(self.ul) // 2, 2)
        #self.ud = np.array(self.ud, dtype='int32')#.reshape(len(self.ud) // 2, 2)

        #if self.k is None or (isinstance(self.k, list) and len(self.k) == 0):
            #self.k = None
        #else:
            #self.k = np.array(self.k)

        #if self.z is None or (isinstance(self.z, list) and len(self.z) == 0):
            #self.z = None
        #else:
            #self.z = np.array(self.z)

        #if self.s is None or (isinstance(self.s, list) and len(self.s) == 0):
            #self.s = None
        #else:
            #self.s = np.array(self.s)

    #def slice_card_by_element_id(self, ids: np.ndarray) -> GENEL:
        #assert self.n > 0, self.n
        #assert len(self.element_id) > 0, self.element_id
        #i = self.index(ids)
        #cls_obj = self.slice_card_by_index(i)
        #assert cls_obj.n > 0, cls_obj
        #return cls_obj

    def __apply_slice__(self, element: GENEL, i: np.ndarray) -> None:
        #print('set_id', self.set_id)
        #print('is_skin', self.is_skin)
        #print('num_ids', self.num_ids)
        #print('ids', self.ids)
        #print('i = ', i)
        #self.set_id = np.zeros(ncards, dtype='int32')
        #self.is_skin = np.zeros(ncards, dtype='bool')
        #self.num_ids = np.zeros(ncards, dtype='int32')
        #self.ids = np.array([], dtype='int32')

        element.element_id = self.element_id[i]

        element.k = hslice_by_idim(i, self.idim_k, self.k)
        element.s = hslice_by_idim(i, self.idim_s, self.s)
        element.z = hslice_by_idim(i, self.idim_z, self.z)

        element.ul = hslice_by_idim(i, self.idim_ul, self.ul)
        element.ud = hslice_by_idim(i, self.idim_ud, self.ud)

        element.nk = self.nk[i]
        element.ns = self.ns[i]
        element.nz = self.nz[i]
        element.nul = self.nul[i]
        element.nud = self.nud[i]

        element.n = len(self.element_id)
        #print('--------------------------------------')
        #print(self)
        assert element.n > 0, element.element_id

    def add(self, eid: int, pid: int, nids: list[int],
            theta_mcid: int|float=0.0, zoffset: float=0.,
            tflag: int=0, T1=None, T2=None, T3=None,
            comment: str=''):
        self.cards.append(((eid, pid, nids,
                            theta_mcid, zoffset,
                            tflag, T1, T2, T3,
                            comment)))
        self.n += 1

    def add_card(self, card: BDFCard, comment=''):
        eid = integer(card, 1, 'eid')
        card_fields = card.fields()
        ucard_fields = [field.upper() if field is not None else None
                        for field in card_fields]
        nfields = card.nfields

        ul = []
        ud = []
        n_ud = ucard_fields.count('UD')

        nk = ucard_fields.count('K')
        nz = ucard_fields.count('Z')
        ns = ucard_fields.count('S')
        assert n_ud in [0, 1], 'n_UD=%s fields=%s' % (n_ud, card_fields)

        assert nk in [0, 1], 'n_K=%s fields=%s' % (nz, card_fields)
        assert nz in [0, 1], 'n_Z=%s fields=%s' % (nk, card_fields)
        assert ns in [0, 1], 'n_S=%s fields=%s' % (ns, card_fields)

        i_ul = 3
        _ul_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_ul)
        for i, _ul in enumerate(_ul_fields):
            uli = integer(card, i + i_ul, 'UL_%d' % (i + 1))
            ul.append(uli)

        if n_ud:
            i_ud = ucard_fields.index('UD')
            if i_ud < nfields and ucard_fields[i_ud] == 'UD':
                assert ucard_fields[i_ud] == 'UD', fields
                _ud_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_ud+2)
                for i, _ud in enumerate(_ud_fields):
                    udi = integer(card, i + i_ud+2, 'UD_%d' % (i + 1))
                    ud.append(udi)

        k = []
        z = []
        s = []
        if nk:
            k = []
            ik = ucard_fields.index('K')
            assert ucard_fields[ik] == 'K', card_fields
            _k_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, ik+1)
            for i, _k in enumerate(_k_fields):
                ki = double(card, i + ik+1, 'K_%d' % (i + 1))
                k.append(ki)
            unused_nblanks = _get_genel_offset(nk)
            #kz = k

        if nz:
            assert len(k) == 0, k
            iz = card_fields.index('Z')
            assert card_fields[iz] == 'Z', card_fields
            _z_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, iz+1)
            for i, _z in enumerate(_z_fields):
                zi = double(card, i + iz+1, 'Z_%d' % (i + 1))
                z.append(zi)
            unused_nblanks = _get_genel_offset(nz)
            #kz = z

        if ns:
            s = []
            i_s = ucard_fields.index('S')
            assert ucard_fields[i_s] == 'S', card_fields
            _s_fields, unused_istop = _read_genel_fields_until_char_blank(ucard_fields, i_s+1)
            for i, _s in enumerate(_s_fields):
                si = double(card, i + i_s+1, 'S_%d' % (i + 1))
                s.append(si)
            unused_nblanks = _get_genel_offset(ns)

        #---------------------------------
        ul = np.array(ul).reshape(len(ul) // 2, 2)
        ud = np.array(ud).reshape(len(ud) // 2, 2)

        #return GENEL(eid, ul, ud, k, z, s, comment=comment)
        self.cards.append((eid, ul, ud, k, z, s, comment))
        self.n += 1

    @Element.parse_cards_check
    def parse_cards(self):
        ncards = len(self.cards)
        idtype = self.model.idtype
        element_id = np.zeros(ncards, dtype=idtype)
        nul = np.zeros(ncards, dtype='int32')
        nud = np.zeros(ncards, dtype='int32')
        nk = np.zeros(ncards, dtype='int32')
        ns = np.zeros(ncards, dtype='int32')
        nz = np.zeros(ncards, dtype='int32')

        ul_list = []
        ud_list = []
        k_list = []
        s_list = []
        z_list = []
        for icard, card in enumerate(self.cards):
            (eid, uli, udi, ki, zi, si, comment) = card
            ns[icard] = len(si)
            nud[icard] = len(udi)
            nul[icard] = len(uli)
            nz[icard] = len(zi)
            nk[icard] = len(ki)
            #print(len(ki), len(si), len(zi))

            element_id[icard] = eid
            k_list.append(ki)
            s_list.append(si)
            z_list.append(zi)
            ul_list.append(uli)
            ud_list.append(udi)
        k = hstack_empty(k_list, dtype_default='float64')
        z = hstack_empty(z_list, dtype_default='float64')
        s = hstack_empty(s_list, dtype_default='float64')
        ul = vstack_empty(ul_list, default_shape=(0, 2), dtype_default='int32')
        ud = vstack_empty(ud_list, default_shape=(0, 2), dtype_default='int32')
        self._save(element_id, ns, nud, nul, nz, nk, k, z, s, ul, ud)
        self.cards = []

    def _save(self, element_id, ns, nud, nul, nz, nk, k, z, s, ul, ud) -> None:
        if len(self.element_id) != 0:
           asdf
        self.element_id = element_id
        self.ns = ns
        self.nud = nud
        self.nul = nul
        self.nz = nz
        self.nk = nk
        self.k = k
        self.z = z
        self.s = s
        self.ul = ul
        self.ud = ud

    #def cross_reference(self, model: BDF) -> None:
        #"""
        #Cross links the card so referenced cards can be extracted directly

        #Parameters
        #----------
        #model : BDF()
            #the BDF object
        #"""
        #msg = ', which is required by GENEL eid=%s' % self.eid
        #self.ul_nodes_ref = model.Nodes(self.ul[:, 0], msg=msg)
        #if len(self.ud):
            #self.ud_nodes_ref = model.Nodes(self.ud[:, 0], msg=msg)

    #def safe_cross_reference(self, model: BDF, xref_errors):
        #"""
        #Cross links the card so referenced cards can be extracted directly

        #Parameters
        #----------
        #model : BDF()
            #the BDF object
        #"""
        ##msg = ', which is required by GENEL eid=%s' % self.eid
        #self.cross_reference(model)
        ##self.ga_ref = model.Node(self.ga, msg=msg)
        ##self.gb_ref = model.Node(self.gb, msg=msg)

    #def uncross_reference(self) -> None:
        #"""Removes cross-reference links"""
        #self.ul[:, 0] = self.ul_nodes
        #self.ud[:, 0] = self.ud_nodes
        #self.ul_nodes_ref = None
        #self.ud_nodes_ref = None

    #def center_of_mass(self):
        #return 0.0

    #@property
    #def nodes(self):
        #return self.node_ids

    #@property
    #def node_ids(self):
        #nodes = self.ul[:, 0].tolist()
        #if len(self.ud):
            #nodes += self.ud[:, 0].tolist()
        #return nodes

    #@property
    #def ul_nodes(self):
        #"""gets the {UL} nodes"""
        #if self.ul_nodes_ref is None:
            #return self.ul[:, 0]
        #nodes = _node_ids(self, nodes=self.ul_nodes_ref, allow_empty_nodes=False)
        #return nodes

    #@property
    #def ud_nodes(self):
        #"""gets the {UD} nodes"""
        #if self.ud_nodes_ref is None:
            #return self.ud[:, 0]
        #nodes = _node_ids(self, nodes=self.ud_nodes_ref, allow_empty_nodes=False)
        #return nodes

    @property
    def idim_k(self) -> np.ndarray:
        return make_idim(self.n, self.nk)
    @property
    def idim_s(self) -> np.ndarray:
        return make_idim(self.n, self.ns)
    @property
    def idim_z(self) -> np.ndarray:
        return make_idim(self.n, self.nz)

    @property
    def idim_ul(self) -> np.ndarray:
        return make_idim(self.n, self.nul)
    @property
    def idim_ud(self) -> np.ndarray:
        return make_idim(self.n, self.nud)

    @property
    def max_id(self) -> int:
        return max(self.element_id.max(), self.ul.max(), self.ud.max())

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike,
                   size: int=8, is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        element_ids = array_str(self.element_id, size=size)

        for eid, (ik0, ik1), (is0, is1), (iz0, iz1), \
            (iul0, iul1), (iud0, iud1), \
            nk, ns, nz in zip(element_ids, \
                              self.idim_k, self.idim_s, self.idim_z,
                              self.idim_ul, self.idim_ud,
                              self.nk, self.ns, self.nz):

            k = self.k[ik0:ik1]
            s = self.s[is0:is1]
            z = self.z[iz0:iz1]
            ud = self.ud[iud0:iud1, :]
            ul = self.ul[iul0:iul1, :]
            list_fields = ['GENEL', eid]

            ## we add 2 to represent the GENEL,eid fields
            n_ul = ul.shape[0] * 2 + 2
            ul_nones = _get_genel_offset(n_ul) * [None]

            ## same as UL for UD
            n_ud = self.ud.shape[0] * 2 + 2
            ud_nones = _get_genel_offset(n_ud) * [None]

            # we call this kz to simplify our life
            kz_char = 'K' if nk is not None else 'Z'
            kz = k if nk is not None else z

            # K/Z has a +1 instead of +2 because there is no blank after K/Z
            n_kz = len(kz) + 1
            kz_nones = _get_genel_offset(n_kz) * [None]

            ud_line = []
            ud_nodes = ud[:, 0]
            ud_dofs = ud[:, 1]
            ud_nodes_dofs = []
            for ud_node, ud_dof in zip(ud_nodes, ud_dofs):
                ud_nodes_dofs.extend([ud_node, ud_dof])
            ud_line = ['UD', None] + ud_nodes_dofs + ud_nones

            #print('s = %r' % self.s, self.s is not None)
            s_line = []
            if ns:
                s_line = ['S'] + s.tolist()

            ul_nodes = ul[:, 0]
            ul_dofs = ul[:, 1]
            ul_nodes_dofs = []
            for ul_node, ul_dof in zip(ul_nodes, ul_dofs):
                ul_nodes_dofs.extend([ul_node, ul_dof])
            ul_line = ul_nodes_dofs + ul_nones

            list_fields = ['GENEL', eid, None] + (
                ul_line +
                ud_line +
                [kz_char] + kz.tolist() + kz_nones +
                s_line
            )
            bdf_file.write(print_card(list_fields))
        return


def _get_genel_offset(n_ul):
    """we add to to represent the GENEL,eid fields"""
    n_ul_leftover = n_ul % 8
    return 8 - n_ul_leftover


def _read_genel_fields_until_char_blank(card_fields, istart):
    """somewhat loose parser helper function for GENEL"""
    new_fields = []
    i = 0
    for i, field in enumerate(card_fields[istart:]):
        if field is None:
            break
        if field.upper() in ['UD', 'K', 'S', 'Z']:
            break
        new_fields.append(field)
    return new_fields, istart+i

def hstack_empty(mylist: list[np.ndarray],
                 dtype_default: str='int32') -> np.ndarray:
    nlist = len(mylist)
    if nlist == 0:
        return np.array([], dtype=dtype_default)
    elif nlist == 1:
        return np.array(mylist[0], dtype=dtype_default)
    else:
        return np.hstack(mylist)

def vstack_empty(mylist: list[np.ndarray],
                 default_shape=None,
                 dtype_default: str='int32') -> np.ndarray:
    nlist = len(mylist)
    if nlist == 0:
        return np.zeros(default_shape, dtype=dtype_default)
    elif nlist == 1:
        return np.array(mylist[0], dtype=dtype_default)
    else:
        return np.vstack(mylist)
