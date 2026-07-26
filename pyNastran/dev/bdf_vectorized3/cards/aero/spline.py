from __future__ import annotations
from itertools import zip_longest
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.field_writer_8 import set_blank_if_default # print_float_8,
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double, double_or_blank,
    string_or_blank, parse_components)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, hslice_by_idim,
    parse_check, save_ifile_comment)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float,
    array_float_nan, get_print_card_size)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.femutils.utils import hstack_lists


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


SPLINE1_MSG = """
+---------+-------+-------+------+------+------+----+------+-------+
|    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
+=========+=======+=======+======+======+======+====+======+=======+
| SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
+---------+-------+-------+------+------+------+----+------+-------+
|         | NELEM | MELEM |      |      |      |    |      |       |
+---------+-------+-------+------+------+------+----+------+-------+
| SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
+---------+-------+-------+------+------+------+----+------+-------+""".strip()

class SPLINE1(VectorizedBaseCard):
    """
    Surface Spline Methods
    Defines a surface spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

    +---------+-------+-------+------+------+------+----+------+-------+
    |    1    |   2   |    3  |   4  |   5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+======+======+======+====+======+=======+
    | SPLINE1 | EID   | CAERO | BOX1 | BOX2 | SETG | DZ | METH | USAGE |
    +---------+-------+-------+------+------+------+----+------+-------+
    |         | NELEM | MELEM |      |      |      |    |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+
    | SPLINE1 |   3   |  111  | 115  | 122  |  14  | 0. |      |       |
    +---------+-------+-------+------+------+------+----+------+-------+

    """
    _id_name = 'spline_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spline_id = np.array([], dtype='int32')
        self.caero_id = np.array([], dtype='int32')

    def add(self, eid: int, caero: int, box1: int, box2: int, setg: int,
            dz: float=0., method: str='IPS',
            usage: str='BOTH',
            nelements: int=10, melements: int=10,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a SPLINE1, which defines a surface spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements : int; default=10
            The number of FE elements along the local spline x-axis if
            using the FPS option
        melements : int; default=10
            The number of FE elements along the local spline y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        card = (eid, caero, box1, box2, setg, dz, method, usage,
                nelements, melements, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPLINE1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', default=0.0)
        method = string_or_blank(card, 7, 'method', default='IPS')
        usage = string_or_blank(card, 8, 'usage', default='BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', default=10)
        melements = integer_or_blank(card, 10, 'melements', default=10)
        assert len(card) <= 11, f'len(SPLINE1 card) = {len(card):d}\ncard={card}'
        #return SPLINE1(eid, caero, box1, box2, setg, dz, method, usage,
                       #nelements, melements, comment=comment)
        card = (eid, caero, box1, box2, setg, dz, method, usage,
                nelements, melements, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        #igroup = np.zeros(ncards, dtype='int32')
        set_id = np.zeros(ncards, dtype='int32')
        box_id = np.zeros((ncards, 2), dtype='int32')
        dz = np.zeros(ncards, dtype='float64')

        usage = np.zeros(ncards, dtype='|U4')
        method = np.zeros(ncards, dtype='|U4')
        nelement = np.zeros(ncards, dtype='int32')
        melement = np.zeros(ncards, dtype='int32')
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, caero, box1, box2, setg, dzi, methodi, usagei,
                nelementi, melementi, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            spline_id[icard] = eid
            caero_id[icard] = caero
            box_id[icard, :] = [box1, box2]
            set_id[icard] = setg
            dz[icard] = dzi
            method[icard] = methodi
            usage[icard] = usagei
            nelement[icard] = nelementi
            melement[icard] = melementi
        self._save(spline_id, caero_id, box_id, set_id, dz,
                   method, usage,
                   nelement, melement,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, set_id, dz,
              method, usage,
              nelement, melement,
              ifile=None, comment=None):
        ncards = len(spline_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.spline_id):
            ifile = np.stack([self.ifile, ifile])
            afd
        save_ifile_comment(self, ifile, comment)
        self.spline_id = spline_id
        self.caero_id = caero_id
        self.box_id = box_id
        self.set_id = set_id
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelement = nelement
        self.melement = melement
        self.n = len(spline_id)

    def __apply_slice__(self, elem: SPLINE1, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.box_id = self.box_id[i, :]
        elem.set_id = self.set_id[i]
        elem.dz = self.dz[i]
        elem.method = self.method[i]
        elem.usage = self.usage[i]
        elem.nelement = self.nelement[i]
        elem.melement = self.melement[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(model.caero1.element_id, ucaero_ids),
        )

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.set_id.max(),
                   self.box_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxes = array_str(self.box_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        nelements = array_default_int(self.nelement, default=0, size=size)
        melements = array_default_int(self.melement, default=0, size=size)
        dzs = array_default_float(self.dz, default=0., size=size, is_double=False)
        for eid, caero, (box1, box2), setg, dz, \
            method, usage, nelement, melement in zip(
                spline_ids, caero_ids, boxes, set_ids, dzs,
                self.method, self.usage, nelements, melements):
            #dz = set_blank_if_default(self.dz, 0.)
            #method = set_blank_if_default(self.method, 'IPS')
            #usage = set_blank_if_default(self.usage, 'BOTH')
            #nelements = set_blank_if_default(self.nelements, 10)
            #melements = set_blank_if_default(self.melements, 10)

            list_fields = ['SPLINE1', eid, caero, box1, box2,
                           setg, dz, method, usage, nelement, melement]
            bdf_file.write(print_card(list_fields))
        return


class SPLINE2(VectorizedBaseCard):
    """
    Linear Spline
    Defines a beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by regular
    arrays of aerodynamic points.

    +---------+------+-------+-------+-------+------+----+------+-----+
    |    1    |   2  |   3   |   4   |   5   |  6   |  7 |   8  |  9  |
    +=========+======+=======+=======+=======+======+====+======+=====+
    | SPLINE2 | EID  | CAERO |  ID1  |  ID2  | SETG | DZ | DTOR | CID |
    +---------+------+-------+-------+-------+------+----+------+-----+
    |         | DTHX | DTHY  | None  | USAGE |      |    |      |     |
    +---------+------+-------+-------+-------+------+----+------+-----+
    | SPLINE2 |   5  |   8   |  12   | 24    | 60   | 0. | 1.0  |  3  |
    +---------+------+-------+-------+-------+------+----+------+-----+
    |         |  1.  |       |       |       |      |    |      |     |
    +---------+------+-------+-------+-------+------+----+------+-----+

    """
    _id_name = 'spline_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spline_id = np.array([], dtype='int32')
        self.caero_id = np.array([], dtype='int32')

    #def add_spline1(self, eid: int, caero: int, box1: int, box2: int, setg: int,
                    #dz: float=0., method: str='IPS',
                    #usage: str='BOTH', nelements: int=10,
                    #melements: int=10, comment: str='') -> int:
        #"""
        #Creates a SPLINE1, which defines a surface spline.

        #Parameters
        #----------
        #eid : int
            #spline id
        #caero : int
            #CAEROx id that defines the plane of the spline
        #box1 / box2 : int
            #First/last box id that is used by the spline
        #setg : int
            #SETx id that defines the list of GRID points that are used
            #by the surface spline
        #dz : float; default=0.0
            #linear attachment flexibility
            #dz = 0.; spline passes through all grid points
        #method : str; default=IPS
            #method for spline fit
            #valid_methods = {IPS, TPS, FPS}
            #IPS : Harder-Desmarais Infinite Plate Spline
            #TPS : Thin Plate Spline
            #FPS : Finite Plate Spline
        #usage : str; default=BOTH
            #Spline usage flag to determine whether this spline applies
            #to the force transformation, displacement transformation, or
            #both
            #valid_usage = {FORCE, DISP, BOTH}
        #nelements : int; default=10
            #The number of FE elements along the local spline x-axis if
            #using the FPS option
        #melements : int; default=10
            #The number of FE elements along the local spline y-axis if
            #using the FPS option
        #comment : str; default=''
            #a comment for the card

        #"""
        #card = (eid, caero, box1, box2, setg, dz, method, usage,
                #nelements, melements, ifile, comment)
        #self.cards.append(card)
        #self.n += 1

    def add(self, eid: int, caero: int,
            box1: int, box2: int, setg: int,
            dz: float=0.0, dtor: float=1.0,
            cid: int=0,
            dthx: float=0.0, dthy: float=0.0,
            usage: str='BOTH',
            ifile: int=0, comment: str='') -> int:
        """
        Creates a SPLINE2 card, which defines a beam spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box/body id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the beam spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        dtor : float; default=1.0
            Torsional flexibility ratio (EI/GJ).
            Use 1.0 for bodies (CAERO2).
        cid : int; default=0
            Rectangular coordinate system for which the y-axis defines the
            axis of the spline. Not used for bodies, CAERO2
        dthx : float; default=None
            Rotational attachment flexibility.
            DTHX : Used for rotation about the spline's x-axis (in-plane
                   bending rotations).  It is not used for bodies (CAERO2).
            DTHY : Used for rotation about the spline's y-axis (torsion).
                   It is used for slope of bodies.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        assert dthx is not None, dthx
        assert dthy is not None, dthy
        card = (eid, caero, box1, box2, setg, dz, dtor, cid,
                dthx, dthy, usage, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPLINE2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box1 = integer(card, 3, 'box1')
        box2 = integer(card, 4, 'box2')
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', default=0.0)
        dtor = double_or_blank(card, 7, 'dtor', default=1.0)
        cid = integer_or_blank(card, 8, 'cid', default=0)
        dthx = double_or_blank(card, 9, 'dthx', default=0.)
        dthy = double_or_blank(card, 10, 'dthy', default=0.)

        usage = string_or_blank(card, 12, 'usage', default='BOTH')
        assert len(card) <= 13, f'len(SPLINE2) card = {len(card):d}\ncard={card}'
        #return SPLINE2(eid, caero, box1, box2, setg, dz, dtor, cid,
                       #dthx, dthy, usage, comment=comment)
        card = (eid, caero, box1, box2, setg, dz, dtor, cid,
                dthx, dthy, usage, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        coord_id = np.zeros(ncards, dtype='int32')
        set_id = np.zeros(ncards, dtype='int32')
        box_id = np.zeros((ncards, 2), dtype='int32')
        dz = np.zeros(ncards, dtype='float64')
        dtor = np.zeros(ncards, dtype='float64')
        dthx = np.zeros(ncards, dtype='float64')
        dthy = np.zeros(ncards, dtype='float64')
        usage = np.zeros(ncards, dtype='|U4')
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, caero, box1, box2, setg, dzi, dtori, cid,
                dthxi, dthyi, usagei, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            spline_id[icard] = eid
            caero_id[icard] = caero
            coord_id[icard] = cid
            box_id[icard, :] = [box1, box2]
            set_id[icard] = setg
            dz[icard] = dzi
            dtor[icard] = dtori
            dthx[icard] = dthxi
            dthy[icard] = dthyi
            usage[icard] = usagei

        self._save(spline_id, caero_id, box_id, set_id, dz,
                   coord_id, dtor, dthx, dthy, usage,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, set_id, dz,
                   coord_id, dtor, dthx, dthy, usage,
              ifile=None, comment=None):
        ncards = len(spline_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.spline_id):
            ifile = np.stack([self.ifile, ifile])
            spline_id = np.hstack([self.spline_id, spline_id])
            caero_id = np.hstack([self.caero_id, caero_id])
            box_id = np.vstack([self.box_id, box_id])
            set_id = np.hstack([self.set_id, set_id])
            coord_id = np.hstack([self.coord_id, coord_id])
            dz = np.hstack([self.dz, dz])
            dtor = np.hstack([self.dtor, dtor])
            usage = np.hstack([self.usage, usage])
            dthx = np.hstack([self.dthx, dthx])
            dthy = np.hstack([self.dthy, dthy])
        save_ifile_comment(self, ifile, comment)
        self.spline_id = spline_id
        self.caero_id = caero_id
        self.box_id = box_id
        self.set_id = set_id
        self.coord_id = coord_id
        self.dz = dz
        self.dtor = dtor
        self.usage = usage
        self.dthx = dthx
        self.dthy = dthy
        self.n = len(spline_id)

    def __apply_slice__(self, elem: SPLINE2, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.ifile = self.ifile[i]
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.box_id = self.box_id[i, :]
        elem.set_id = self.set_id[i]
        elem.coord_id = self.coord_id[i]

        elem.dthx = self.dthx[i]
        elem.dthy = self.dthy[i]
        elem.dz = self.dz[i]
        elem.dtor = self.dtor[i]
        #elem.dtor = self.dtor[i]
        elem.usage = self.usage[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        all_caero_ids = hstack_lists([model.caero1.element_id, model.caero2.element_id],
                                     unique_sort=True)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(all_caero_ids, ucaero_ids),
        )

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.set_id.max(),
                   self.box_id.max(), self.coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxes = array_str(self.box_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        dzs = array_default_float(self.dz, default=0.0, size=size, is_double=False)
        dtors = array_default_float(self.dtor, default=1.0, size=size, is_double=False)
        dthxs = array_default_float(self.dthx, default=0.0, size=size, is_double=False)
        dthys = array_default_float(self.dthy, default=0.0, size=size, is_double=False)
        for eid, caero, (box1, box2), setg, dz, dtor, \
            cid, usage, dthx, dthy in zip(
                spline_ids, caero_ids, boxes, set_ids, dzs, dtors,
                coord_ids, self.usage, dthxs, dthys):

            list_fields = ['SPLINE2', eid, caero, box1, box2,
                           setg, dz, dtor, cid, dthx, dthy,
                           None, usage]
            bdf_file.write(print_card(list_fields))
        return


class SPLINE3(VectorizedBaseCard):
    """
    Defines a constraint equation for aeroelastic problems.
    Useful for control surface constraints.

    +---------+------+-------+-------+------+----+----+-----+-------+
    |    1    |  2   |   3   |   4   |  5   |  6 |  7 |  8  |   9   |
    +=========+======+=======+=======+======+====+====+=====+=======+
    | SPLINE3 | EID  | CAERO | BOXID | COMP | G1 | C1 | A1  | USAGE |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  G2  |  C2   |  A2   |      | G3 | C3 | A2  |       |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  G4  |  C4   |  A4   | etc. |    |    |     |       |
    +---------+------+-------+-------+------+----+----+-----+-------+
    | SPLINE3 | 7000 |  107  |  109  |  6   | 5  | 3  | 1.0 | BOTH  |
    +---------+------+-------+-------+------+----+----+-----+-------+
    |         |  43  |   5   | -1.0  |      |    |    |     |       |
    +---------+------+-------+-------+------+----+----+-----+-------+

    """
    _id_name = 'spline_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spline_id = np.array([], dtype='int32')
        self.caero_id = np.array([], dtype='int32')

        self.components = np.array([], dtype='int32')
        self.box_id = np.array([], dtype='int32')

        self.nnode = np.array([], dtype='int32')
        self.usage = np.array([], dtype='|U4')
        self.nodes = np.array([], dtype='int32')
        self.displacement_components = np.array([], dtype='int32')
        self.coeffs = np.array([], dtype='float64')

    def add(self, eid: int, caero: int, box_id: int,
            components: int,
            nodes: list[int],
            displacement_components: list[int],
            coeffs: list[float],
            usage: str='BOTH',
            ifile: int=0, comment: str='') -> int:
        """
        Creates a SPLINE3 card, which is useful for control surface
        constraints.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box_id : int
           Identification number of the aerodynamic box number.
        components : int
           The component of motion to be interpolated.
           3, 5          (CAERO1)
           2, 3, 5, 6    (CAERO2)
           3             (CAERO3)
           3, 5, 6       (CAERO4)
           3, 5, 6       (CAERO5)
           1, 2, 3, 5, 6 (3D Geometry)
           2-lateral displacement
           3-transverse displacement
           5-pitch angle
           6-relative control angle for CAERO4/5; yaw angle for CAERO2

        nodes : list[int]
           Grid point identification number of the independent grid point.
        displacement_components : list[int]
           Component numbers in the displacement coordinate system.
           1-6 (GRIDs)
           0 (SPOINTs)
        coeffs : list[float]
           Coefficient of the constraint relationship.
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        comment : str; default=''
            a comment for the card

        """
        nnode = len(nodes)
        assert nnode == len(displacement_components)
        assert nnode == len(coeffs)
        card = (eid, caero, box_id, components,
                nodes, displacement_components,
                coeffs, usage, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPLINE3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        box_id = integer(card, 3, 'box_id')
        components = integer(card, 4, 'comp')
        node = integer(card, 5, 'G1')
        displacement_component = integer(card, 6, 'C1')
        coeff = double(card, 7, 'A1')
        usage = string_or_blank(card, 8, 'usage', default='BOTH')

        nfields = len(card) - 1
        nrows = nfields // 8
        if nfields % 8:
            nrows += 1

        nodes = [node]
        coeffs = [coeff]
        displacement_components = [displacement_component]
        i = 2
        for irow in range(1, nrows):
            #print('G%i' % i)
            j = 1 + irow * 8
            node = integer(card, j, 'G%i' % i)
            displacement_component = integer(card, j + 1, 'C%i' % i)
            coeff = double(card, j + 2, 'A%i' % i)
            nodes.append(node)
            coeffs.append(coeff)
            displacement_components.append(displacement_component)
            i += 1
            if card.field(j + 4) or card.field(j + 5) or card.field(j + 6):
                node = integer(card, j + 4, 'G%i' % i)
                displacement_component = parse_components(card, j + 5, 'C%i' % i)
                coeff = double(card, j + 6, 'A%i' % i)
                nodes.append(node)
                coeffs.append(coeff)
                displacement_components.append(int(displacement_component))
                i += 1
        #spline = SPLINE3(eid, caero, box_id, components,
                         #nodes, displacement_components, coeffs, usage=usage,
                         #comment=comment)
        card = (eid, caero, box_id, components,
                nodes, displacement_components,
                coeffs, usage, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        components = np.zeros(ncards, dtype='int32')
        box_id = np.zeros(ncards, dtype='int32')

        nnode = np.zeros(ncards, dtype='int32')
        #ndisp = np.zeros(ncards, dtype='int32')
        #dz = np.zeros(ncards, dtype='float64')
        usage = np.zeros(ncards, dtype='|U4')
        comment = {}

        nodes = []
        displacement_components = []
        coeffs = []
        for icard, card in enumerate(self.cards):
            (eid, caero, box_idi, componentsi,
                nodesi, displacement_componentsi,
                coeffsi, usagei, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            spline_id[icard] = eid
            caero_id[icard] = caero
            box_id[icard] = box_idi
            #set_id[icard] = setg
            #dz[icard] = dzi
            components[icard] = componentsi
            usage[icard] = usagei
            nnode[icard] = len(nodesi)
            nodes.extend(nodesi)
            displacement_components.extend(displacement_componentsi)
            coeffs.extend(coeffsi)

        nodes = np.array(nodes, dtype='int32')
        displacement_components = np.array(displacement_components, dtype='int32')
        coeffs = np.array(coeffs, dtype='float64')
        self._save(spline_id, caero_id, box_id, components,
                nnode, nodes, displacement_components,
                coeffs, usage, ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, components,
                nnode, nodes, displacement_components, coeffs, usage,
              ifile=None, comment=None):
        ncards = len(spline_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.spline_id):
            ifile = np.stack([self.ifile, ifile])
            spline_id = np.hstack([self.spline_id, spline_id])
            caero_id = np.hstack([self.caero_id, caero_id])
            box_id = np.hstack([self.box_id, box_id])
            components = np.hstack([self.components, components])
            nnode = np.hstack([self.nnode, nnode])
            nodes = np.hstack([self.nodes, nodes])
            displacement_components = np.hstack([self.displacement_components, displacement_components])
            coeffs = np.hstack([self.coeffs, coeffs])
            usage = np.hstack([self.usage, usage])

        assert len(nodes) > 0
        assert len(displacement_components) > 0
        assert len(coeffs) > 0
        save_ifile_comment(self, ifile, comment)
        self.spline_id = spline_id
        self.caero_id = caero_id
        self.box_id = box_id
        self.components = components

        self.nnode = nnode
        self.nodes = nodes
        self.displacement_components = displacement_components
        self.coeffs = coeffs
        self.usage = usage
        self.n = len(spline_id)

    #def sort(self) -> None:
        #ueid = np.unique(self.spline_id)
        #if np.array_equal(ueid, self.spline_id):
            #return
        #i = np.argsort(self.spline_id)
        #self.__apply_slice__(self, i)

    def __apply_slice__(self, elem: SPLINE3, i: np.ndarray) -> None:
        self.write()
        elem.n = len(i)
        elem.ifile = self.ifile[i]
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.box_id = self.box_id[i]
        elem.usage = self.usage[i]
        elem.components = self.components[i]

        inode = self.inode
        elem.displacement_components = hslice_by_idim(i, inode, self.displacement_components)
        elem.coeffs = hslice_by_idim(i, inode, self.coeffs)
        elem.nodes = hslice_by_idim(i, inode, self.nodes)
        elem.nnode = self.nnode[i]
        self.write()

    def validate(self):
        msg = ''
        cmin = self.components.min()
        cmax = self.components.max()
        components_allowed = {0, 1, 2, 3, 4, 5, 6}
        if cmax not in components_allowed:
            msg += f'components.min()={cmin} must be [0, 1, 2, 3, 4, 5, 6]\n'
        if cmin not in components_allowed:
            msg += f'components.max()={cmax} must be [0, 1, 2, 3, 4, 5, 6]\n'

        #if not len(self.nodes) == len(self.displacement_components):
            #msg += 'nnodes=%s ndisplacement_components=%s must be equal\n' % (
                #len(self.nodes), len(self.displacement_components))
        #if not len(self.nodes) == len(self.coeffs):
            #msg += 'nnodes=%s ncoeffs=%s must be equal\n' % (
                #len(self.nodes), len(self.coeffs))

        for i, disp_component  in enumerate(self.displacement_components):
            if disp_component not in components_allowed:
                if not isinstance(disp_component, integer_types):
                    msg += (
                        f'i={i} displacement_component={disp_component!r} must be an integer '
                        f'[0, 1, 2, 3, 4, 5, 6]; type={type(disp_component)}\n')
                else:
                    msg += f'i={i} displacement_component={disp_component} must be [0, 1, 2, 3, 4, 5, 6]\n'
        for usage in np.unique(self.usage):
            if usage not in ['FORCE', 'DISP', 'BOTH']:
                msg += f'usage={usage} must be in [FORCE, DISP, BOTH]\n'

        if msg:
            msg += str(self)
            raise RuntimeError(msg)

        #for node in self.nodes:
            #assert isinstance(node, integer_types), self.nodes
        #for displacement_component in self.displacement_components:
            #assert isinstance(displacement_component, integer_types), self.displacement_components
        #for coeff in self.coeffs:
            #assert isinstance(coeff, float), self.coeffs

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(model.caero1.element_id, ucaero_ids),
        )

    @property
    def inode(self) -> np.ndarray:
        inode = make_idim(self.n, self.nnode)
        return inode

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.box_id.max())

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.box_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxes = array_str(self.box_id, size=size)
        components = array_str(self.components, size=size)
        nodes = array_str(self.nodes, size=size)
        disp = array_str(self.displacement_components, size=size)

        #spline_id, caero_id, box_id, components,
        #nnodes, nodes, displacement_components, coeffs
        for eid, caero, box_id, component, (inode0, inode1), usage in zip(
                spline_ids, caero_ids, boxes, components,
                self.inode, self.usage):
            usages = set_blank_if_default(usage, 'BOTH')
            nodesi = nodes[inode0:inode1]
            dispi = disp[inode0:inode1]
            coeffsi = self.coeffs[inode0:inode1]

            list_fields = [
                'SPLINE3', eid, caero, box_id, component,
                nodesi[0], dispi[0], coeffsi[0], usages]
            for nid, dispii, coeff in zip(nodes[1:], dispi[1:], coeffsi[1:]):
                list_fields += [nid, dispii, coeff, None]

            bdf_file.write(print_card(list_fields))
        return


class SPLINE4(VectorizedBaseCard):
    """
    Surface Spline Methods
    Defines a curved surface spline for interpolating motion and/or forces for
    aeroelastic problems on general aerodynamic geometries using either the
    Infinite Plate, Thin Plate or Finite Plate splining method.

    +---------+-------+-------+--------+-----+------+----+------+-------+
    |    1    |   2   |   3   |    4   |  5  |   6  |  7 |   8  |   9   |
    +=========+=======+=======+========+=====+======+====+======+=======+
    | SPLINE4 |  EID  | CAERO | AELIST |     | SETG | DZ | METH | USAGE |
    +---------+-------+-------+--------+-----+------+----+------+-------+
    |         | NELEM | MELEM |        |     |      |    |      |       |
    +---------+-------+-------+--------+-----+------+----+------+-------+
    | SPLINE4 |   3   | 111   |   115  |     |  14  | 0. | IPS  |       |
    +---------+-------+-------+--------+-----+------+----+------+-------+

    """
    _id_name = 'spline_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spline_id = np.array([], dtype='int32')
        self.caero_id = np.array([], dtype='int32')

    def add(self, eid: int, caero: int, aelist: int, setg: int,
            dz: float=0., method: str='IPS', usage: str='BOTH',
            nelements: int=10, melements: int=10,
            ftype: Optional[str]='WF2',
            rcore: Optional[float]=None,
            ifile: int=0, comment: str='') -> int:
        """
        Creates a SPLINE4 card, which defines a curved Infinite Plate,
        Thin Plate, or Finite Plate Spline.

        Parameters
        ----------
        eid : int
            spline id
        caero : int
            CAEROx id that defines the plane of the spline
        box1 / box2 : int
            First/last box id that is used by the spline
        setg : int
            SETx id that defines the list of GRID points that are used
            by the surface spline
        dz : float; default=0.0
            linear attachment flexibility
            dz = 0.; spline passes through all grid points
        method : str; default=IPS
            method for spline fit
            valid_methods = {IPS, TPS, FPS}
            IPS : Harder-Desmarais Infinite Plate Spline
            TPS : Thin Plate Spline
            FPS : Finite Plate Spline
        usage : str; default=BOTH
            Spline usage flag to determine whether this spline applies
            to the force transformation, displacement transformation, or
            both
            valid_usage = {FORCE, DISP, BOTH}
        nelements / melements : int; default=10
            The number of FE elements along the local spline x/y-axis if
            using the FPS option
        comment : str; default=''
            a comment for the card

        """
        card = (eid, caero, aelist, setg, dz, method, usage,
                nelements, melements, ftype, rcore,
                ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setg')
        dz = double_or_blank(card, 6, 'dz', default=0.0)
        method = string_or_blank(card, 7, 'method', default='IPS')
        usage = string_or_blank(card, 8, 'usage', default='BOTH')
        nelements = integer_or_blank(card, 9, 'nelements', default=10)
        melements = integer_or_blank(card, 10, 'melements', default=10)
        ftype = string_or_blank(card, 11, 'ftype', default='WF2')
        rcore = double_or_blank(card, 12, 'rcore', default=np.nan)
        assert len(card) <= 13, f'len(SPLINE4) card = {len(card):d}\ncard={card}'
        #return SPLINE4(eid, caero, aelist, setg, dz, method, usage,
                       #nelements, melements, ftype=ftype, rcore=rcore,
                       #comment=comment)
        card = (eid, caero, aelist, setg, dz, method, usage,
                nelements, melements, ftype, rcore,
                ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        aelist_id = np.zeros(ncards, dtype='int32')
        set_id = np.zeros(ncards, dtype='int32')
        dz = np.zeros(ncards, dtype='float64')

        usage = np.zeros(ncards, dtype='|U4')
        method = np.zeros(ncards, dtype='|U4')
        nelement = np.zeros(ncards, dtype='int32')
        melement = np.zeros(ncards, dtype='int32')
        ftype = np.zeros(ncards, dtype='|U4')
        rcore = np.zeros(ncards, dtype='float64')
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, caero, aelist, setg, dzi, methodi, usagei,
                nelementi, melementi, ftypei, rcorei,
                ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            spline_id[icard] = eid
            caero_id[icard] = caero
            aelist_id[icard] = aelist
            set_id[icard] = setg
            dz[icard] = dzi
            method[icard] = methodi
            usage[icard] = usagei
            nelement[icard] = nelementi
            melement[icard] = melementi
            ftype[icard] = ftypei
            rcore[icard] = rcorei
        self._save(spline_id, caero_id, aelist_id, set_id, dz,
                   method, usage,
                   nelement, melement, ftype, rcore,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, aelist_id, set_id, dz,
              method, usage,
              nelement, melement, ftype, rcore,
              ifile=None, comment=None):
        ncards = len(spline_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.spline_id):
            ifile = np.stack([self.ifile, ifile])
            afd
        save_ifile_comment(self, ifile, comment)
        self.spline_id = spline_id
        self.caero_id = caero_id
        self.aelist_id = aelist_id
        self.set_id = set_id
        self.dz = dz
        self.method = method
        self.usage = usage
        self.nelement = nelement
        self.melement = melement
        self.ftype = ftype
        self.rcore = rcore
        self.n = len(spline_id)

    def __apply_slice__(self, elem: SPLINE4, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.ifile = self.ifile[i]
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.aelist_id = self.aelist_id[i]
        elem.set_id = self.set_id[i]
        elem.dz = self.dz[i]
        elem.method = self.method[i]
        elem.usage = self.usage[i]
        elem.nelement = self.nelement[i]
        elem.melement = self.melement[i]
        elem.ftype = self.ftype[i]
        elem.rcore = self.rcore[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(model.caero1.element_id, ucaero_ids),
        )

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.set_id.max(),
                   self.aelist_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        aelist_ids = array_str(self.aelist_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        nelements = array_default_int(self.nelement, default=0, size=size)
        melements = array_default_int(self.melement, default=0, size=size)
        dzs = array_float(self.dz, size=size)
        rcores = array_float_nan(self.rcore, size=size, is_double=False)
        for eid, caero, aelist, setg, dz, \
            method, usage, nelement, melement, ftype, rcore in zip(
                spline_ids, caero_ids, aelist_ids, set_ids, dzs,
                self.method, self.usage, nelements, melements, self.ftype, rcores):

            list_fields = ['SPLINE4', eid, caero, aelist, None,
                           setg, dz, method, usage, nelement, melement,
                           ftype, rcore]
            bdf_file.write(print_card(list_fields))
        return


class SPLINE5(VectorizedBaseCard):
    """
    Linear Spline
    Defines a 1D beam spline for interpolating motion and/or forces for
    aeroelastic problems on aerodynamic geometries defined by irregular arrays
    of aerodynamic points. The interpolating beam supports axial rotation and
    bending in the yz-plane.

    +=========+======+=======+========+=======+======+====+=======+=======+
    |    1    |  2   |    3  |    4   |   5   |   6  |  7 |   8   |   9   |
    +=========+======+=======+========+=======+======+====+=======+=======+
    | SPLINE5 | EID  | CAERO | AELIST |       | SETG | DZ | DTOR  |  CID  |
    +---------+------+-------+--------+-------+------+----+-------+-------+
    |         | DTHX | DTHY  |        | USAGE | METH |    | FTYPE | RCORE |
    +---------+------+-------+--------+-------+------+----+-------+-------+

    METH, FTYPE, RCORE are in 2012+ (not MSC.2005r2 or NX.10)

    """
    _id_name = 'spline_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.spline_id = np.array([], dtype='int32')
        self.caero_id = np.array([], dtype='int32')

        self.aelist_id = np.array([], dtype='int32')
        self.set_id = np.array([], dtype='int32')
        self.dz = np.array([], dtype='float64')
        self.dtor = np.array([], dtype='float64')
        self.coord_id = np.array([], dtype='int32')

        self.usage = np.array([], dtype='|U4')
        self.method = np.array([], dtype='|U4')
        self.ftype = np.array([], dtype='|U4')
        self.rcore = np.array([], dtype='float64')

    def add(self, eid: int, caero: int, aelist: int, setg: int, thx, thy,
            dz: float=0.0, dtor: float=1.0, cid: int=0,
            usage: str='BOTH', method: str='BEAM',
            ftype: str='WF2', rcore=None,
            ifile: int=0, comment: str='') -> int:
        """Creates a SPLINE5 card"""
        card = (eid, caero, aelist, setg, thx, thy, dz, dtor, cid,
                usage, method, ftype, rcore, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, ifile: int, comment: str='') -> int:
        """
        Adds a SPLINE5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        caero = integer(card, 2, 'caero')
        aelist = integer(card, 3, 'aelist')
        # None
        setg = integer(card, 5, 'setq')
        dz = double_or_blank(card, 6, 'dz', default=0.0)
        dtor = double_or_blank(card, 7, 'dtor', default=1.0)
        cid = integer_or_blank(card, 8, 'cid', default=0)
        thx = double(card, 9, 'thx')
        thy = double(card, 10, 'thy')
        usage = string_or_blank(card, 12, 'usage', default='BOTH')
        # per nast/tpl/fmondsp.dat, METH can be a double(0.0) ???
        method = string_or_blank(card, 13, 'meth', default='BEAM')
        ftype = string_or_blank(card, 15, 'ftype', default='WF2')
        rcore = double_or_blank(card, 16, 'rcore')

        usage = string_or_blank(card, 12, 'usage', default='BOTH')
        assert len(card) <= 16, 'len(SPLINE5 card) = %i\n%s' % (len(card), card)
        #return SPLINE5(eid, caero, aelist, setg, thx, thy, dz=dz, dtor=dtor, cid=cid,
                       #usage=usage, method=method, ftype=ftype, rcore=rcore, comment=comment)
        card = (eid, caero, aelist, setg, thx, thy, dz, dtor, cid,
                usage, method, ftype, rcore, ifile, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        ifile = np.zeros(ncards, dtype='int32')
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        aelist_id = np.zeros(ncards, dtype='int32')
        set_id = np.zeros(ncards, dtype='int32')
        thx = np.zeros(ncards, dtype='float64')
        thy = np.zeros(ncards, dtype='float64')
        dz = np.zeros(ncards, dtype='float64')
        dtor = np.zeros(ncards, dtype='float64')
        coord_id = np.zeros(ncards, dtype='int32')

        usage = np.zeros(ncards, dtype='|U4')
        method = np.zeros(ncards, dtype='|U4')
        ftype = np.zeros(ncards, dtype='|U4')
        rcore = np.zeros(ncards, dtype='float64')
        comment = {}
        for icard, card in enumerate(self.cards):
            (eid, caero, aelist, setg, thxi, thyi, dzi, dtori, cidi,
             usagei, methodi, ftypei, rcorei, ifilei, commenti) = card
            ifile[icard] = ifilei
            if commenti:
                comment[eid] = commenti
            spline_id[icard] = eid
            caero_id[icard] = caero
            aelist_id[icard] = aelist
            set_id[icard] = setg
            thx[icard] = thxi
            thy[icard] = thyi
            dz[icard] = dzi
            dtor[icard] = dtori
            coord_id[icard] = cidi
            method[icard] = methodi
            usage[icard] = usagei
            ftype[icard] = ftypei
            rcore[icard] = rcorei
        self._save(spline_id, caero_id, aelist_id, set_id,
                   thx, thy, dz, dtor,
                   coord_id, method, usage, ftype, rcore,
                   ifile=ifile, comment=comment)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, aelist_id, set_id,
              thx, thy, dz, dtor,
              coord_id, method, usage, ftype, rcord,
              ifile=None, comment=None):
        ncards = len(spline_id)
        if ifile is None:
            ifile = np.zeros(ncards, dtype='int32')
        if len(self.spline_id):
            ifile = np.stack([self.ifile, ifile])
            afd
        save_ifile_comment(self, ifile, comment)
        self.spline_id = spline_id
        self.caero_id = caero_id
        self.aelist_id = aelist_id
        self.set_id = set_id
        self.thx = thx
        self.thy = thy
        self.dz = dz
        self.dtor = dtor
        self.coord_id = coord_id
        self.method = method
        self.usage = usage
        self.ftype = ftype
        self.rcord = rcord
        self.n = len(spline_id)

    def __apply_slice__(self, elem: SPLINE5, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.aelist_id = self.aelist_id[i]
        elem.set_id = self.set_id[i]
        elem.thx = self.thx[i]
        elem.thy = self.thy[i]
        elem.dz = self.dz[i]
        elem.dtor = self.dtor[i]
        elem.coord_id = self.coord_id[i]
        elem.method = self.method[i]
        elem.usage = self.usage[i]
        elem.ftype = self.ftype[i]
        elem.rcord = self.rcord[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(model.caero1.element_id, ucaero_ids),
        )

    @property
    def max_id(self) -> int:
        return max(self.spline_id.max(), self.caero_id.max(), self.set_id.max(),
                   self.coord_id.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card, size = get_print_card_size(size, self.max_id)

        spline_ids = array_str(self.spline_id, size=size)
        caero_ids = array_str(self.caero_id, size=size)
        aelists = array_str(self.aelist_id, size=size)
        cids = array_str(self.coord_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        thxs = array_float(self.thx, size=size, is_double=False)
        thys = array_float(self.thy, size=size, is_double=False)
        dzs = array_default_float(self.dz, default=0., size=size, is_double=False)
        for eid, caero, aelist, setg, thx, thy, dz, dtor, \
            cid, method, usage, ftype, rcore in zip_longest(
                spline_ids, caero_ids, aelists, set_ids, thxs, thys, dzs, self.dtor,
                cids, self.method, self.usage, self.ftype, self.rcore):

            #dz = set_blank_if_default(dz, 0.)
            usage = set_blank_if_default(usage, 'BOTH')

            list_fields = ['SPLINE5', eid, caero, aelist, None,
                           setg, dz, dtor, cid, thx, thy,
                           None, usage, method, None, ftype, rcore]
            bdf_file.write(print_card(list_fields))
        return
