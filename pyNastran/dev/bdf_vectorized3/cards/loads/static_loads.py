from __future__ import annotations
from collections import defaultdict
from typing import Union, TYPE_CHECKING
import numpy as np

from pyNastran.bdf.cards.base_card import expand_thru_by
from pyNastran.bdf.field_writer_8 import set_string8_blank_if_default, print_card_8, print_float_8 # , print_field_8
from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, double, string,
    integer_or_blank, double_or_blank, integer_or_string,
    components_or_blank, fields)
from pyNastran.bdf.cards.collpase_card import collapse_thru_by
from pyNastran.bdf.bdf_interface.assign_type_force import force_integer
from pyNastran.utils.numpy_utils import (
    integer_types, float_types,   # integer_float_types,
)

from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.dev.bdf_vectorized3.cards.coord import transform_spherical_to_rectangular
from pyNastran.dev.bdf_vectorized3.cards.base_card import VectorizedBaseCard, hslice_by_idim, make_idim # , searchsorted_filter
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int

if TYPE_CHECKING:
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from .dynamic_loads import LOADSET


class Load(VectorizedBaseCard):
    _id_name = 'load_id'
    def __init__(self, model: BDF):
        super().__init__(model)
        self.load_id = np.array([], dtype='int32')

    def slice_card_by_load_id(self, load_id: np.ndarray) -> Load0:
        assert len(self.load_id) > 0, self
        card_class = self.__class__
        card = card_class(self.model)
        if isinstance(load_id, integer_types):
            load_id_set = {load_id}
        else:
            load_id_set = set(load_id.tolist())

        index = []
        for i, load_idi in enumerate(self.load_id):
            if load_idi in load_id_set:
                index.append(i)
        assert len(index) > 0, f'no {card.type}s found; {self} load_id={load_id} load_ids={self.load_id}'
        return self.slice_card_by_index(index)


class Load0(Load):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.node_id = np.array([], dtype='int32')
        self.coord_id = np.array([], dtype='int32')
        self.mag = np.array([], dtype='float64')
        self.xyz = np.zeros((0, 3), dtype='float64')

    def add(self, sid: int, node: int, mag: float, xyz: np.ndarray,
            cid: int=0, comment: str='') -> FORCE:
        """
        Creates a FORCE/MOMENT card

        Parameters
        ----------
        sid : int
            load id
        node : int
            the node to apply the load to
        mag : float
            the load's magnitude
        xyz : (3, ) float ndarray
            the load direction in the cid frame
        cid : int; default=0
            the coordinate system for the load
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, node, cid, mag, xyz, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        cid = integer_or_blank(card, 3, 'cid', default=0)
        mag = double(card, 4, 'mag')
        xyz = [double_or_blank(card, 5, 'X1', default=0.0),
               double_or_blank(card, 6, 'X2', default=0.0),
               double_or_blank(card, 7, 'X3', default=0.0)]
        assert len(card) <= 8, 'len(%s card) = %d\ncard=%s' % (self.type, len(card), card)
        self.cards.append((sid, node, cid, mag, xyz, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        load_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')
        mag = np.zeros(ncards, dtype='float64')
        xyz = np.zeros((ncards, 3), dtype='float64')
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, node, cid, magi, xyzi, comment) = card
            load_id[icard] = sid
            node_id[icard] = node
            coord_id[icard] = cid
            mag[icard] = magi
            xyz[icard, :] = xyzi
        self._save(load_id, node_id, coord_id, mag, xyz)
        assert len(load_id) == self.n
        self.cards = []

    def _save(self, load_id, node_id, coord_id, mag, xyz):
        if len(self.load_id) != 0:
            load_id = np.hstack([self.load_id, load_id])
            node_id = np.hstack([self.node_id, node_id])
            coord_id = np.hstack([self.coord_id, coord_id])
            mag = np.hstack([self.mag, mag])
            xyz = np.vstack([self.xyz, xyz])
        nloads = len(load_id)
        self.load_id = load_id
        self.node_id = node_id
        self.coord_id = coord_id
        self.mag = mag
        self.xyz = xyz
        self.n = nloads

    def geom_check(self, missing: dict[str, np.ndarray]):
        nid = self.model.grid.node_id
        cid = self.model.coord.coord_id

        #load_nodes = self.node_id
        #load_cid = self.coord_id
        #assert base_nodes is not None
        #print(self.base_nodes)
        geom_check(self,
                   missing,
                   node=(nid, self.node_id), filter_node0=False,
                   coord=(cid, self.coord_id))

    def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        if len(self.load_id) == 0:
            return ''
        lines = []
        card_class = self.type
        if size == 8:
            for sid, nid, cid, mag, xyz in zip(self.load_id, self.node_id, self.coord_id, self.mag, self.xyz):
                cids = set_string8_blank_if_default(cid, 0)
                msg = '%-8s%8d%8d%8s%8s%8s%8s%8s\n' % (
                    card_class, sid, nid,
                    cids, print_float_8(mag), print_float_8(xyz[0]),
                    print_float_8(xyz[1]), print_float_8(xyz[2]))
                lines.append(msg)
        else:
            raise RuntimeError(size)
        return ''.join(lines)

class Load1(Load):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.node_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 2), dtype='float64')
        self.mag = np.array([], dtype='float64')
        #self.xyz = np.zeros((0, 3), dtype='float64')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        assert len(card) == 6, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        self.cards.append((sid, node, mag, [g1, g2], comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        load_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        mag = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 2), dtype='int32')
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            sid, node, magi, g12, comment = card
            load_id[icard] = sid
            node_id[icard] = node
            mag[icard] = magi
            nodes[icard, :] = g12
        self._save(load_id, node_id, mag, nodes)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, node_id, mag, nodes):
        if len(self.load_id) != 0:
            raise NotImplementedError()
        nloads = len(load_id)
        self.load_id = load_id
        self.node_id = node_id
        self.mag = mag
        self.nodes = nodes
        self.n = nloads

    def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        if len(self.load_id) == 0:
            return ''
        lines = []
        card_class = self.type
        if size == 8:
            for sid, nid, mag, nodes in zip(self.load_id, self.node_id, self.mag, self.nodes):
                msg = '%-8s%8d%8d%8s%8s%8s\n' % (
                    card_class, sid, nid, print_float_8(mag), nodes[0], nodes[1])
                lines.append(msg)
        else:
            raise RuntimeError(size)
        return ''.join(lines)

class Load2(Load):
    def __init__(self, model: BDF):
        super().__init__(model)
        self.node_id = np.array([], dtype='int32')
        self.nodes = np.zeros((0, 4), dtype='float64')
        self.mag = np.array([], dtype='float64')
        #self.xyz = np.zeros((0, 3), dtype='float64')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        node = integer(card, 2, 'node')
        mag = double(card, 3, 'mag')
        g1 = integer(card, 4, 'g1')
        g2 = integer(card, 5, 'g2')
        g3 = integer(card, 6, 'g3')
        g4 = integer(card, 7, 'g4')
        assert len(card) == 8, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        self.cards.append((sid, node, mag, [g1, g2, g3, g4], comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        ncards = len(self.cards)
        if ncards == 0:
            return
        load_id = np.zeros(ncards, dtype='int32')
        node_id = np.zeros(ncards, dtype='int32')
        mag = np.zeros(ncards, dtype='float64')
        nodes = np.zeros((ncards, 4), dtype='int32')
        assert ncards > 0, ncards
        for icard, card in enumerate(self.cards):
            (sid, node, magi, nodesi, comment) = card
            load_id[icard] = sid
            node_id[icard] = node
            mag[icard] = magi
            nodes[icard, :] = nodesi
        self._save(load_id, node_id, mag, nodes)
        assert len(self.load_id) == self.n
        self.cards = []

    def _save(self, load_id, node_id, mag, nodes):
        if len(self.load_id) != 0:
            raise NotImplementedError()
        nloads = len(load_id)
        self.load_id = load_id
        self.node_id = node_id
        self.mag = mag
        self.nodes = nodes
        self.n = nloads

    def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        if len(self.load_id) == 0:
            return ''
        lines = []
        card_class = self.type
        load_ids = array_str(self.load_id, size=size)
        node_id_ = array_str(self.node_id, size=size)
        node_ids_ = array_default_int(self.nodes, default=0, size=size)
        if size == 8:
            for sid, nid, mag, nodes in zip(load_ids, node_id_, self.mag, node_ids_):
                msg = '%-8s%8s%8s%8s%8s%8s%8s%8s\n' % (
                    card_class, sid, nid, print_float_8(mag), nodes[0], nodes[2], nodes[1], nodes[3])
                lines.append(msg)
        else:
            raise RuntimeError(size)
        return ''.join(lines)

class FORCE(Load0):
    def slice_card_by_index(self, i: np.ndarray) -> FORCE:
        load = FORCE(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: FORCE, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.coord_id = self.coord_id[i]
        load.mag = self.mag[i]
        load.xyz = self.xyz[i, :]

    @property
    def scaled_vector(self) -> np.ndarray:
        return self.mag * self.xyz

    def sum_forces_moments(self) -> np.ndarray:
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()
        nid = grid.node_id

        nloads = len(self.load_id)
        force_moment = np.zeros((nloads, 6), dtype='float64')
        force = force_moment[:, :3]
        moment = force_moment[:, 3:]

        ucoords = np.unique(self.coord_id)
        for cid in ucoords:
            icoord = np.where(self.coord_id == cid)[0]
            force0 = self.mag[icoord, np.newaxis] * self.xyz[icoord, :]
            moment0 = np.cross(force0, self.xyz[icoord, :])
            #print('cid =', cid)
            #print('xyz =', self.xyz[icoord, :])
            #print('force0', force0)
            #print('moment0', moment0)
            if cid == 0:
                force[icoord, :] = force0
                moment[icoord, :] = moment0
            else:
                #print('else...')
                coord = self.model.coord
                coord_card = coord.slice_card_by_coord_id(cid)
                beta = coord_card.xyz_to_global_transform[cid]

                coord_type = coord_card.coord_type
                if coord_type == 'R':
                    # TODO: I'm pretty sure this is right...
                    force[icoord, :] = force0 @ beta
                    moment[icoord, :] = moment0 @ beta
                elif coord_type == 'S':
                    # TODO: I'm pretty sure this is right...
                    force_r = np.array([transform_spherical_to_rectangular(force_) for force_ in force0])
                    moment_r = np.array([transform_spherical_to_rectangular(moment_) for moment_ in moment0])
                    force[icoord, :] = force_r @ beta
                    moment[icoord, :] = moment_r @ beta
                else:
                    raise NotImplementedError(f'coord={cid} not supported\n{coord_card.write()}')
        return force_moment

class MOMENT(Load0):
    def slice_card_by_index(self, i: np.ndarray) -> MOMENT:
        load = MOMENT(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: MOMENT, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.coord_id = self.coord_id[i]
        load.mag = self.mag[i]
        load.xyz = self.xyz[i, :]

    @property
    def scaled_vector(self) -> np.ndarray:
        return self.mag * self.xyz

    def sum_forces_moments(self):
        grid = self.model.grid
        #xyz_cid0 = grid.xyz_cid0()
        #nid = grid.node_id

        nloads = len(self.load_id)
        ucoords = np.unique(self.coord_id)
        moment = self.mag[:, np.newaxis] * self.xyz
        #if self.coord_id.max() == 0 and self.coord_id.min() == 0:
            #moment = local_moment
        #else:

        coord =  self.model.coord
        for ucid in ucoords:
            if ucid == 0:
                continue
            icid = np.where(ucid == self.coord_id)[0]
            local_moment = moment[icid, :]
            global_moment = coord.transform_force_local_to_global(local_moment)
            moment[icid, :] = global_moment
            #coords = coord.slice_card_by_coord_id(ucoords)
            #raise NotImplementedError(f'the following coordinate systems are not supported\n{coords.write()}')

        force_moment = np.zeros((nloads, 6), dtype='float64')
        force_moment[:, 3:] = moment
        return force_moment


class FORCE1(Load1):
    """
    Defines a static concentrated force at a grid point by specification of a
    magnitude and two grid points that determine the direction.

    +--------+-----+----+-------+----+----+
    |   1    |  2  | 3  |   4   | 5  | 6  |
    +========+=====+====+=======+====+====+
    | FORCE1 | SID | G  |   F   | G1 | G2 |
    +--------+-----+----+-------+----+----+
    | FORCE1 |  6  | 13 | -2.93 | 16 | 13 |
    +--------+-----+----+-------+----+----+

    """
    def slice_card_by_index(self, i: np.ndarray) -> FORCE1:
        load = FORCE1(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: FORCE1, i: np.ndarray) -> None:  # ignore[override]
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.mag = self.mag[i]
        load.nodes = self.nodes[i, :]

    def sum_forces_moments(self):
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        nloads = len(self.load_id)
        iapplied_nid = np.searchsorted(grid.node_id, self.node_id)
        inid = np.searchsorted(grid.node_id, self.nodes)
        in1 = inid[:, 0]
        in2 = inid[:, 1]

        xyz = xyz_cid0[iapplied_nid, :]
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        nxyz = xyz2 - xyz1
        dist = np.linalg.norm(nxyz, axis=1)

        assert len(dist) == nloads
        assert dist.min() > 0, dist
        nxyz /= dist[:, np.newaxis]

        force = self.mag[:, np.newaxis] * nxyz
        moment = np.cross(xyz, force)

        force_moment = np.hstack([force, moment])
        assert force_moment.shape == (nloads, 6)
        return force_moment

class MOMENT1(Load1):
    def slice_card_by_index(self, i: np.ndarray) -> MOMENT1:
        load = MOMENT1(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: MOMENT1, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.mag = self.mag[i]
        load.nodes = self.nodes[i, :]

    def sum_forces_moments(self):
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        nloads = len(self.load_id)
        iapplied_nid = np.searchsorted(grid.node_id, self.node_id)
        inid = np.searchsorted(grid.node_id, self.nodes)
        in1 = inid[:, 0]
        in2 = inid[:, 1]

        xyz = xyz_cid0[iapplied_nid, :]
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        nxyz = xyz2 - xyz1
        dist = np.linalg.norm(nxyz, axis=1)

        assert len(dist) == nloads
        assert dist.min() > 0, dist
        nxyz /= dist[:, np.newaxis]

        moment = self.mag[:, np.newaxis] * nxyz
        force_moment = np.hstack([np.zeros((nloads, 3), dtype='float64'), moment])
        assert force_moment.shape == (nloads, 6)
        return force_moment


class FORCE2(Load2):
    def slice_card_by_index(self, i: np.ndarray) -> FORCE2:
        load = FORCE2(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: FORCE2, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.mag = self.mag[i]
        load.nodes = self.nodes[i, :]

    def sum_forces_moments(self):
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        nloads = len(self.load_id)
        iapplied_nid = np.searchsorted(grid.node_id, self.node_id)
        inid = np.searchsorted(grid.node_id, self.nodes)
        in1 = inid[:, 0]
        in2 = inid[:, 1]
        in3 = inid[:, 2]
        in4 = inid[:, 3]
        assert in4.min() > 0, in4

        xyz = xyz_cid0[iapplied_nid, :]
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        xyz3 = xyz_cid0[in3, :]
        xyz4 = xyz_cid0[in4, :]

        #v21 = xyz2 - xyz1
        #v2 = xyz4 - xyz3
        #xyz = cross(v21, v2)
        nxyz = np.cross(xyz2 - xyz1, xyz4 - xyz3)
        dist = np.linalg.norm(nxyz, axis=1)

        assert len(dist) == nloads
        assert dist.min() > 0, dist
        nxyz /= dist[:, np.newaxis]

        force = self.mag[:, np.newaxis] * nxyz
        moment = np.cross(xyz, force)

        force_moment = np.hstack([force, moment])
        assert force_moment.shape == (nloads, 6)
        return force_moment

class MOMENT2(Load2):
    def slice_card_by_index(self, i: np.ndarray) -> MOMENT2:
        load = MOMENT2(self.model)
        self.__apply_slice__(load, i)
        return load

    def __apply_slice__(self, load: MOMENT2, i: np.ndarray) -> None:
        load.n = len(i)
        load.load_id = self.load_id[i]
        load.node_id = self.node_id[i]
        load.mag = self.mag[i]
        load.nodes = self.nodes[i, :]

    def sum_forces_moments(self):
        grid = self.model.grid
        xyz_cid0 = grid.xyz_cid0()

        nloads = len(self.load_id)
        iapplied_nid = np.searchsorted(grid.node_id, self.node_id)
        inid = np.searchsorted(grid.node_id, self.nodes)
        in1 = inid[:, 0]
        in2 = inid[:, 1]
        in3 = inid[:, 2]
        in4 = inid[:, 3]
        assert in4.min() > 0, in4

        xyz = xyz_cid0[iapplied_nid, :]
        xyz1 = xyz_cid0[in1, :]
        xyz2 = xyz_cid0[in2, :]
        xyz3 = xyz_cid0[in3, :]
        xyz4 = xyz_cid0[in4, :]

        #v21 = xyz2 - xyz1
        #v2 = xyz4 - xyz3
        #xyz = cross(v21, v2)
        nxyz = np.cross(xyz2 - xyz1, xyz4 - xyz3)
        dist = np.linalg.norm(nxyz, axis=1)

        assert len(dist) == nloads
        assert dist.min() > 0, dist
        nxyz /= dist[:, np.newaxis]

        moment = self.mag[:, np.newaxis] * nxyz
        force_moment = np.hstack([np.zeros((nloads, 3), dtype='float64'), moment])
        assert force_moment.shape == (nloads, 6)
        return force_moment


class LOAD(Load):
    """
    +------+-----+------+------+----+-----+----+----+----+
    |   1  |  2  |  3   |  4   | 5  |  6  | 7  | 8  | 9  |
    +======+=====+======+======+====+=====+====+====+====+
    | LOAD | SID |  S   |  S1  | L1 | S2  | L2 | S3 | L3 |
    +------+-----+------+------+----+-----+----+----+----+
    |      | S4  |  L4  | etc. |    |     |    |    |    |
    +------+-----+------+------+----+-----+----+----+----+
    | LOAD | 101 | -0.5 | 1.0  | 3  | 6.2 | 4  |    |    |
    +------+-----+------+------+----+-----+----+----+----+

    """
    def __init__(self, model: BDF):
        super().__init__(model)
        self.nloads = np.array([], dtype='int32')
        self.load_ids = np.array([], dtype='int32')
        self.scale_factors = np.array([], dtype='float64')

    def add(self, sid: int, scale: float,
            scale_factors: list[float],
            load_ids: list[int], comment: str='') -> None:
        assert len(scale_factors) == len(load_ids), f'sid={sid:d} scale_factors={scale_factors} load_ids={load_ids}'
        self.cards.append((sid, scale, scale_factors, load_ids, comment))
        self.n += 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        scale = double(card, 2, 'scale')

        scale_factors = []
        load_ids = []

        # alternating of scale factor & load set ID
        nload_fields = len(card) - 3
        assert nload_fields % 2 == 0, 'card=%s' % card
        for iload in range(nload_fields // 2):
            n = 2 * iload + 3
            scale_factors.append(double(card, n, 'scale_factor'))
            load_ids.append(integer(card, n + 1, 'load_id'))

        assert len(card) > 3, 'len(%s card) = %i\ncard=%s' % (self.type, len(card), card)
        self.cards.append((sid, scale, scale_factors, load_ids, comment))
        self.n += 1
        return self.n

    def parse_cards(self) -> None:
        if self.n == 0:
            return
        nloads = len(self.cards)
        if nloads == 0:
            return
        self.load_id = np.zeros(nloads, dtype='int32')
        self.scale = np.zeros(nloads, dtype='float64')
        self.nloads = np.zeros(nloads, dtype='int32')

        all_load_ids = []
        all_scale_factors = []
        assert nloads > 0, nloads
        for icard, card in enumerate(self.cards):
            (sid, scale, scale_factors, load_ids, comment) = card

            nloads_actual = len(scale_factors)

            self.load_id[icard] = sid
            self.scale[icard] = scale
            self.nloads[icard] = nloads_actual
            all_load_ids.extend(load_ids)
            all_scale_factors.extend(scale_factors)
        self.load_ids = np.array(all_load_ids, dtype='int32')
        self.scale_factors = np.array(all_scale_factors, dtype='float64')
        self.cards = []

    @property
    def iload(self) -> np.ndarray:
        return make_idim(self.n, self.nloads)

    def write(self, size: int=8, is_double: bool=False, write_card_header: bool=False) -> str:
        if len(self.load_id) == 0:
            return ''
        #get_reduced_loads(self, filter_zero_scale_factors=False)
        lines = []

        for sid, scale, iload in zip(self.load_id, self.scale, self.iload):
            iload0, iload1 = iload
            list_fields = ['LOAD', sid, scale]
            scale_factors = self.scale_factors[iload0:iload1]
            load_ids = self.load_ids[iload0:iload1]
            for (scale_factor, load_id) in zip(scale_factors, load_ids):
                list_fields += [scale_factor, load_id]
            #if len(load_ids) != len(scale_factors):
                #msg = 'nload_ids=%s nscale_factors=%s and arent the same\n' % (
                    #len(load_ids), len(scale_factors))
                #msg = 'load_ids=%s\n' % (load_ids)
                #msg += 'scale_factors=%s\n' % (scale_factors)
                #msg += print_card_8(list_fields)
                #msg += str(self.get_stats())
                #raise IndexError(msg)
            lines.append(print_card_8(list_fields))
        #else:
            #raise RuntimeError(size)
        return ''.join(lines)


    #def get_loads_by_load_id(self) -> dict[int, Loads]:
        #return get_loads_by_load_id(self)

    def get_loads_by_load_id(load: Union[LOAD, LOADSET]) -> dict[int, Loads]:
        """"""
        model = load.model
        #uload_ids = np.unique(self.load_ids)
        loads_by_load_id = defaultdict(list)

        #print('all_laods =', model.loads)
        for loadi in model.loads:
            if loadi.type in {'LOAD', 'LSEQ'}:
                continue
            if loadi.n == 0:
                continue
            uload_idsi = np.unique(loadi.load_id)
            #print(f'load.type={loadi.type} {uload_idsi}')
            for uload_id in uload_idsi:
                #print(loadi.load_id)
                #i = np.where(uload_id == loadi.load_id)[0]
                #if len(i) == 0:
                    #print('i =', i)
                    #jj
                    #continue
                #print(f'load.type={loadi.type} {uload_id}; i={i}')
                #loadi = loadi.slice_card_by_index(i)
                loadi2 = loadi.slice_card_by_load_id(uload_id)
                #if loadi.type == 'PLOAD4':
                    #loadi.nvector

                loads_by_load_id[uload_id].append(loadi2)
        return dict(loads_by_load_id)

    def get_reduced_loads(self,
                          remove_missing_loads: bool=False,
                          filter_zero_scale_factors: bool=False,
                          stop_on_failure: bool=True) -> dict[int, Loads]:
        return get_reduced_loads(
            self, remove_missing_loads=remove_missing_loads,
            filter_zero_scale_factors=filter_zero_scale_factors,
            stop_on_failure=stop_on_failure)

    def get_reduced_load_by_load_id(self,
                                    load_id: int,
                                    remove_missing_loads: bool=False,
                                    filter_zero_scale_factors: bool=False,
                                    stop_on_failure: bool=True) -> dict[int, Loads]:
        """
        Parameters
        ----------
        resolve_load_card : bool; default=False
            ???
        remove_missing_loads: bool; default=False
            LOAD cards can reference loads (e.g., GRAV) that don't exist
            Nastran sometimes ignores these loads leading to potentially incorrect results
        filter_zero_scale_factors: bool; default=False
            remove loads that are 0.0
        """
        load = self.slice_card_by_load_id(load_id)
        reduced_loads = get_reduced_loads(load)
        return reduced_loads


def get_reduced_static_load_from_load_id(model: BDF,
                                         load_id: int,
                                         remove_missing_loads: bool=False,
                                         filter_zero_scale_factors: bool=False,
                                         stop_on_failure: bool=True) -> list[StaticLoad]:
    #log = model.log

    load: LOAD = model.load
    reduced_loads = []
    if load.n and load_id in load.load_id:
        load.get_reduced_load_by_load_id(
            load_id,
            remove_missing_loads=remove_missing_loads,
            filter_zero_scale_factors=filter_zero_scale_factors,
            stop_on_failure=stop_on_failure)
        #reduced_loads = load.get_reduced_loads(
            #remove_missing_loads=False,
            #filter_zero_scale_factors=False,
            #stop_on_failure=True)
        raise RuntimeError('aaa')
    else:
        for load in model.loads:
            if load.n == 0:
                continue
            if load.type in {'LOAD', 'LSEQ'}:
                model.log.debug(f'skipping {load.type}')
                continue
            scale_factor = 1.
            loadi = load.slice_card_by_load_id(load_id)
            reduced_loads.append((scale_factor, loadi))
    return reduced_loads

def get_reduced_loads(load: Union[LOAD, LSEQ],
                      remove_missing_loads: bool=False,
                      filter_zero_scale_factors: bool=False,
                      stop_on_failure: bool=True) -> dict[int, Loads]:
    """
    Parameters
    ----------
    resolve_load_card : bool; default=False
        ???
    remove_missing_loads: bool; default=False
        LOAD cards can reference loads (e.g., GRAV) that don't exist
        Nastran sometimes ignores these loads leading to potentially incorrect results
    filter_zero_scale_factors: bool; default=False
        remove loads that are 0.0
    """
    reduced_loads = {}
    if load.n == 0:
        return reduced_loads

    stop_on_failure = True
    loads_by_load_id = load.get_loads_by_load_id()
    log = load.model.log
    for sid, global_scale, iload in zip(load.load_id, load.scale_factors, load.iload):
        reduced_loadsi = []
        iload0, iload1 = iload
        if global_scale == 0. and filter_zero_scale_factors:
            print('continueA')
            continue
        scale_factors = global_scale * load.scale_factors[iload0:iload1]
        load_ids = load.load_ids[iload0:iload1]
        for (scale_factor, load_id) in zip(scale_factors, load_ids):
            if scale_factor == 0. and filter_zero_scale_factors:
                continue
            loads_found = loads_by_load_id[load_id]
            if len(loads_found) == 0:
                msg = f'No referenced loads found for load_id={load_id} on {load.type} load_id={sid}'
                log.error(msg)
                if stop_on_failure:
                    raise RuntimeError(msg)
            reduced_loadsi.append((scale_factor, loads_found))
        reduced_loads[sid] = reduced_loadsi

    # loads that weren't referenced by a LOAD card
    for load_id, loads in loads_by_load_id.items():
        if load_id not in reduced_loads:
            reduced_loads[load_id] = [(1., loads)]
    return reduced_loads
