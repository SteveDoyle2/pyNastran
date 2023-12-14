from __future__ import annotations
from abc import abstractmethod
from itertools import count, zip_longest
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types # , cast_ints
from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default # print_float_8,
#from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, integer_or_string,
    double, double_or_blank,
    string, double_string_or_blank, string_or_blank, # integer_double_or_blank,
    parse_components, # components_or_blank,
    fields, interpret_value)
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default

from pyNastran.bdf.cards.base_card import expand_thru
from pyNastran.bdf.cards.aero.utils import (
    #elements_from_quad,
    points_elements_from_quad_points, # create_ellipse,
    create_axisymmetric_body,
)

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, make_idim, hslice_by_idim,
    get_print_card_8_16,
    parse_element_check,
    parse_property_check,
)
from pyNastran.dev.bdf_vectorized3.cards.write_utils import (
    array_str, array_float,
    array_default_int, array_default_float, array_default_str,
    array_float_nan)
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
from pyNastran.femutils.utils import hstack_lists


if TYPE_CHECKING:
    from pyNastran.nptyping_interface import NDArray3float
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike


class AECOMP(VectorizedBaseCard):
    """
    Defines a component for use in monitor point definition or external splines.

    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |   1    |   2   |    3     |   4   |   5   |   6   |   7   |   8   |   9   |
    +========+=======+==========+=======+=======+=======+=======+=======+=======+
    | AECOMP | NAME  | LISTTYPE | LIST1 | LIST2 | LIST3 | LIST4 | LIST5 | LIST6 |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    |        | LIST7 |   etc.   |       |       |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | WING  |  AELIST  | 1001  | 1002  |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | WING  |   SET1   | 1001  | 1002  |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+
    | AECOMP | WING  |  CAERO1  | 1001  | 2001  |       |       |       |       |
    +--------+-------+----------+-------+-------+-------+-------+-------+-------+

    Attributes
    ----------
    name : str
        The name.
    list_type : str
        {'SET1', 'AELIST', 'CAEROx'}
    lists : list[int]
        list of values of AECOMP lists

    """
    _id_name = 'name'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='|U8')
        self.list_type = np.array([], dtype='|U8')
        self.nlists = np.array([], dtype='int32')
        self.all_lists = np.array([], dtype='int32')

    def add(self, name: str, list_type: list[str], lists: Union[int, list[int]],
            comment: str='') -> int:
        """
        Creates an AECOMP card

        Parameters
        ----------
        name : str
            the name of the component
        list_type : str
            One of CAERO, AELIST or CMPID for aerodynamic components and
            SET1 for structural components. Aerodynamic components are
            defined on the aerodynamic ks-set mesh while the structural
            components are defined on the g-set mesh.
        lists : list[int, int, ...]; int
            The identification number of either SET1, AELIST or CAEROi
            entries that define the set of grid points that comprise
            the component
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((name, list_type, lists, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AECOMP card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        list_type = string(card, 2, 'list_type')
        assert list_type in {'SET1', 'AELIST', 'CAERO'}, list_type # 'CAEROx
        j = 1
        lists = []
        for i in range(3, len(card)):
            list_i = integer(card, i, '%s_%d' % (list_type, j))
            lists.append(list_i)
            j += 1
        self.cards.append((name, list_type, lists, comment))
        self.n += 1
        #return AECOMP(name, list_type, lists, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        list_type = np.zeros(ncards, dtype='|U8')
        nlists = np.zeros(ncards, dtype='int32')

        all_lists = []
        for icard, card in enumerate(self.cards):
            namei, list_typei, listsi, comment = card
            name[icard] = namei
            list_type[icard] = list_typei
            nlists[icard] = len(listsi)
            all_lists.extend(listsi)
        lists = np.array(all_lists, dtype='int32')
        self._save(name, list_type, nlists, lists)
        self.sort()
        self.cards = []

    def _save(self, name, list_type, nlists, lists) -> None:
        assert len(self.name) == 0, self.name
        self.name = name
        self.list_type = list_type
        self.nlists = nlists
        self.lists = lists
        self.n = len(name)

    #def sort(self) -> None:
        #uname = np.unique(self.name)
        #if np.array_equal(uname, self.name):
            #return
        #i = np.argsort(self.name)
        #self.__apply_slice__(self, i)

    def validate(self):
        pass

    def __apply_slice__(self, elem: AECOMP, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.name = self.name[i]
        elem.list_type = self.list_type[i]
        ilist = self.ilist
        elem.lists = hslice_by_idim(i, ilist, self.lists)
        elem.nlists = self.nlists[i]
        elem.n = len(self.name)

    @property
    def ilist(self) -> np.ndarray:
        return make_idim(self.n, self.nlists)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        set1_ids = []
        aelist_ids = []
        caero_ids = []
        for list_type, (ilist0, ilist1) in zip(self.list_type, self.ilist):
            lists = self.lists[ilist0:ilist1]
            if list_type == 'SET1':
                set1_ids.append(lists)
            elif list_type == 'AELIST':
                aelist_ids.append(lists)
            elif list_type == 'CAERO':
                caero_ids.append(lists)
            else:
                raise NotImplementedError(list_type)

        if set1_ids:
            set1_ids = np.unique(np.hstack(set1_ids))
        if aelist_ids:
            aelist_ids = np.unique(np.hstack(aelist_ids))
        if caero_ids:
            caero_ids = np.unique(np.hstack(caero_ids))

        all_caero_ids = model.caero1.element_id
        geom_check(
            self,
            missing,
            set1=(model.set1.set_id, set1_ids),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(all_caero_ids, caero_ids),
        )

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.name) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        lists_ = array_str(self.lists, size=size).tolist()

        for name, list_type, (ilist0, ilist1) in zip(self.name, self.list_type, self.ilist):
            lists = lists_[ilist0:ilist1]
            list_fields = ['AECOMP', name, list_type] + lists
            bdf_file.write(print_card(list_fields))
        return


class AECOMPL(VectorizedBaseCard):
    """
    Makes a "AECOMP" that is a combination of other AECOMPs or AECOMPLs.

    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    |    1    |    2   |    3   |    4   |    5    |    6   |    7   |    8   |    9   |
    +=========+========+========+========+=========+========+========+========+========+
    | AECOMPL |  NAME  | LABEL1 | LABEL2 | LABEL3  | LABEL4 | LABEL5 | LABEL6 | LABEL7 |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    |         | LABEL8 |  etc.  |        |         |        |        |        |        |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    | AECOMPL |  HORIZ |  STAB  |  ELEV  | BALANCE |        |        |        |        |
    +---------+--------+--------+--------+---------+--------+--------+--------+--------+
    """
    _id_name = 'name'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='|U8')
        self.nlabels = np.array([], dtype='int32')
        self.labels = np.array([], dtype='|U8')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AECOMPL card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        name = string(card, 1, 'name')
        labels = []
        j = 1
        for i in range(2, len(card)):
            label = string(card, i, 'label_%d' % j)
            labels.append(label)
            j += 1
        self.cards.append((name, labels, comment))
        self.n += 1
        #return AECOMPL(name, labels, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        nlabels = np.zeros(ncards, dtype='int32')
        #all_labels = np.array([], dtype='|U8')

        all_labels = []
        for icard, card in enumerate(self.cards):
            namei, labelsi, comment = card
            name[icard] = namei
            nlabels[icard] = len(labelsi)
            all_labels.extend(labelsi)
        labels = np.array(all_labels, dtype='|U8')
        self._save(name, nlabels, labels)
        self.sort()
        self.cards = []

    def _save(self, name, nlabels, labels) -> None:
        assert len(self.name) == 0, self.name
        self.name = name
        self.nlabels = nlabels
        self.labels = labels
        self.n = len(name)

    def __apply_slice__(self, elem: AECOMPL, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.name = self.name[i]
        ilabel = self.ilabel
        elem.labels = hslice_by_idim(i, ilabel, self.labels)
        elem.nlabels = self.nlabels[i]

    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabels)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #model = self.model
        for (ilabel0, ilabel1) in self.ilabel:
            labels = self.labels[ilabel0:ilabel1]
        #geom_check(
            #missing,
            ##aecomp=(model.aecomp.aelist_id, aelist_ids),
        #)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.name) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        labels_ = array_str(self.labels, size=size).tolist()
        for name, (ilabel0, ilabel1) in zip(self.name, self.ilabel):
            labels = labels_[ilabel0:ilabel1]
            list_fields = ['AECOMPL', name] + labels
            bdf_file.write(print_card(list_fields))
        return


class CAERO1(VectorizedBaseCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO1 | EID | PID | CP | NSPAN | NCHORD |  LSPAN | LCHORD | IGID |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int
        int : PAERO1 ID
    igroup : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int
        int : coordinate system
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    nchord : int
        int > 0 : N chordwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
    lchord : int
        int > 0 : AEFACT reference for non-uniform nchord
        int = 0 : use nchord
    comment : str; default=''
         a comment for the card

    """
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int, igroup: int,
            p1: NDArray3float, x12: float,
            p4: NDArray3float, x43: float,
            cp: int=0,
            nspan: int=0, lspan: int=0,
            nchord: int=0, lchord: int=0, comment: str='') -> int:
        """
        Defines a CAERO1 card, which defines a simplified lifting surface
        (e.g., wing/tail).

        Parameters
        ----------
        eid : int
            element id
        pid : int
            int : PAERO1 ID
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2; (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3; (typically x, tip chord)
        cp : int; default=0
            int : coordinate system
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        nchord : int; default=0
            int > 0 : N chordwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        lchord : int, AEFACT; default=0
            int > 0 : AEFACT reference for non-uniform nchord
            int = 0 : use nchord
        comment : str; default=''
             a comment for the card

        """
        nchord = nchord if nchord is not None else 0
        lchord = lchord if lchord is not None else 0
        nspan = nspan if nspan is not None else 0
        lspan = lspan if lspan is not None else 0
        cp = cp if cp is not None else 0
        card = (eid, pid, igroup, p1, x12, p4, x43, cp,
                nspan, lspan, nchord, lchord, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        nchord = integer_or_blank(card, 5, 'nchord', default=0)
        lspan = integer_or_blank(card, 6, 'lspan', default=0)
        lchord = integer_or_blank(card, 7, 'lchord', default=0)
        igroup = integer(card, 8, 'igid')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', default=0.0),
            double_or_blank(card, 10, 'y1', default=0.0),
            double_or_blank(card, 11, 'z1', default=0.0)])
        x12 = double_or_blank(card, 12, 'x12', default=0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', default=0.0),
            double_or_blank(card, 14, 'y4', default=0.0),
            double_or_blank(card, 15, 'z4', default=0.0)])
        x43 = double_or_blank(card, 16, 'x43', default=0.)

        assert len(card) <= 17, f'len(CAERO1 card) = {len(card):d}\ncard={card}'
        #return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      #cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      #comment=comment)
        card = (eid, pid, igroup, p1, x12, p4, x43, cp,
                nspan, lspan, nchord, lchord, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        igroup = np.zeros(ncards, dtype='int32')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nspan = np.zeros(ncards, dtype='int32')
        nchord = np.zeros(ncards, dtype='int32')
        lspan = np.zeros(ncards, dtype='int32')
        lchord = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            eid, pid, igroupi, p1i, x12i, p4i, x43i, cpi, nspani, lspani, nchordi, lchordi, comment = card
            element_id[icard] = eid
            property_id[icard] = pid
            igroup[icard] = igroupi
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            nspan[icard] = nspani
            nchord[icard] = nchordi
            lspan[icard] = lspani
            lchord[icard] = lchordi
        self._save(element_id, property_id, igroup, p1, p4, x12, x43, cp,
                   nspan, lspan, nchord, lchord)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, igroup, p1, p4, x12, x43, cp,
              nspan, lspan, nchord, lchord):
        if len(self.element_id):
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            igroup = np.hstack([self.igroup, igroup])
            p1 = np.vstack([self.p1, p1])
            p4 = np.vstack([self.p4, p4])
            x12 = np.hstack([self.x12, x12])
            x43 = np.hstack([self.x43, x43])
            cp = np.hstack([self.cp, cp])
            nspan = np.hstack([self.nspan, nspan])
            lspan = np.hstack([self.lspan, lspan])
            nchord = np.hstack([self.nchord, nchord])
            lchord = np.hstack([self.lchord, lchord])
        self.element_id = element_id
        self.property_id = property_id
        self.igroup = igroup
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.nchord = nchord
        self.lchord = lchord
        self.n = len(element_id)

    def add_quad(self, eid: int, pid: int, span: int, chord: int, igroup: int,
                 p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray,
                 cp: int=0, spanwise: str='y', comment: str='') -> None:
        r"""
        ::

          1
          | \
          |   \
          |     \
          |      4
          |      |
          |      |
          2------3

        TODO: CP not handled correctly

        """
        x12 = p2[0] - p1[0]
        x43 = p3[0] - p4[0]
        nspan = 0
        lspan = 0
        nchord = 0
        lchord = 0

        if spanwise.lower() == 'y':
            y41 = p4[1] - p1[1]
            y32 = p3[1] - p2[1]
            dspan = max(y41, y32)
        elif spanwise.lower() == 'z':
            y41 = p4[2] - p1[2]
            y32 = p3[2] - p2[2]
            dspan = max(y41, y32)
        else:
            raise NotImplementedError('spanwise=%r; expected=[y, z]' % spanwise.lower())

        dx = max(x12, x43)
        if isinstance(span, integer_types):
            nspan = span
        #elif isinstance(span, AEFACT):
            #lspan = span.sid
        elif isinstance(span, float):
            nspan = int(np.ceil(dspan / span))
            if nspan <= 0:
                msg = 'y41=%s y32=%s; dspan=%s span=%s nspan=%s; nspan must be greater than 0' % (
                    y41, y32, dspan, span, nspan)
                raise ValueError(msg)
        else:
            raise TypeError(span)

        if isinstance(chord, integer_types):
            nchord = chord
        #elif isinstance(chord, AEFACT):
            #lchord = chord.sid
        elif isinstance(chord, float):
            nchord = int(np.ceil(dx / chord))
            if nchord <= 0:
                msg = 'x12=%s x43=%s; dx=%s chord=%s nchord=%s; nchord must be greater than 0' % (
                    x12, x43, dx, chord, nchord)
                raise ValueError(msg)
        else:
            raise TypeError(chord)
        n = self.add(eid, pid, igroup, p1, x12, p4, x43, cp=cp,
                 nspan=nspan, lspan=lspan,
                 nchord=nchord, lchord=lchord, comment=comment)
        #return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      #cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      #comment=comment)
        return n

    def __apply_slice__(self, elem: CAERO1, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.igroup = self.igroup[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.nspan = self.nspan[i]
        elem.lspan = self.lspan[i]
        elem.nchord = self.nchord[i]
        elem.lchord = self.lchord[i]

    def validate(self):
        msg = ''
        is_failed = False
        #if not isinstance(self.p1, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True
        #if not isinstance(self.p4, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True

        element_id = self.element_id

        ibad = (self.x12 <= 0.)
        if np.any(ibad):
            msg += 'X12 and must be greater than or equal to 0\n'
            msg += f'  element_id={element_id[ibad]}\n   x12={self.x12[ibad]}'
            is_failed = True
        #if self.x43 <= 0.:
            #msg += 'X43=%s and must be greater than or equal to 0\n' % (self.x43)
            #is_failed = True

        ibad = (self.nspan == 0) & (self.lspan == 0)
        if np.any(ibad):
            msg += 'NSPAN or LSPAN must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nspan != 0) & (self.lspan != 0)
        if np.any(ibad):
            msg += 'Either NSPAN or LSPAN must 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nchord == 0) & (self.lchord == 0)
        if np.any(ibad):
            msg += 'NCHORD or LCHORD must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True

        ibad = (self.nchord != 0) & (self.lchord != 0)
        if np.any(ibad):
            msg += 'Either NCHORD or LCHORD must 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True
        if is_failed:
            msg += str(self)
            #msg += CAERO1_MSG
            raise ValueError(msg)

        assert self.p1.shape[1] == 3, 'p1=%s' % self.p1.shape
        assert self.p4.shape[1] == 3, 'p4=%s' % self.p4.shape

        # calculating area; assuming coordinate transformations don't matter
        if 1:
            p1 = self.p1
            p4 = self.p4
            d12 = np.zeros(p1.shape)
            d12[:, 0] = self.x12
            d43 = np.zeros(p1.shape)
            d43[:, 0] = self.x43
            p2 = p1 + d12
            p3 = p4 + d43

            a = p3 - p1
            b = p4 - p2
            area = np.linalg.norm(np.cross(a, b), axis=1)
            assert len(area) == p1.shape[0]
            ibad = (area < 0.0)
            if np.any(ibad):
                msg += 'Either NCHORD or LCHORD must 0\n'
                msg += f'eid={self.element_id[ibad,:]} p1={p1[ibad,:]} p2={p2[ibad,:]} p3={p3[ibad,:]} p4={p4[ibad,:]} area={area[ibad]}\n'
                raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def flip_normal_by_element_id(self, element_id=None):
        """flips the CAERO1 normal vector"""
        i = self.index(element_id)
        self.flip_normal_by_index(i)

    def flip_normal_by_index(self, i: np.ndarray) -> None:
        """flips the CAERO1 normal vector"""
        self.p1[i, :], self.p4[i, :] = self.p4[i, :], self.p1[i, :]
        self.x12[i, :], self.x43[i, :] = self.x43[i, :], self.x12[i, :]

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        #xy = []
        x = []
        y = []
        aefact = self.model.aefact
        for eid, nchord, lchord, nspan, lspan in zip(self.element_id, self.nchord, self.lchord, self.nspan, self.lspan):
            if nchord == 0:
                aefacti = aefact.slice_card_by_aefact_id(lchord)
                nchord = aefacti.nfractions - 1  # points -> boxes
                xi = aefacti.fractions
            else:
                xi = np.linspace(0., 1., nchord + 1)

            if nspan == 0:
                aefacti = aefact.slice_card_by_aefact_id(lspan)
                nspan = aefacti.nfractions - 1  # points -> boxes
                yi = aefacti.fractions
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nchord < 1 or nspan < 1:
                msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                    eid, nchord, nspan, lchord, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    def get_leading_edge_points(self):
        """gets the leading edge points"""
        p1 = self.p1
        p4 = self.p4
        icps = (self.cp != 0)
        if np.any(icps):
            #coord = self.model.coord.slice_card_by_coord_id(self.cp)
            cps = self.cp[icps]
            for icp, cp in zip(icps, cps):
                p1[icp, :] = self.model.coord.transform_local_xyz_to_global(self.p1[icp, :], cp)
                p4[icp, :] = self.model.coord.transform_local_xyz_to_global(self.p4[icp, :], cp)
        return p1, p4

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        p1, p4 = self.get_leading_edge_points()
        zero = np.zeros(self.x12.size, dtype=self.x12.dtype)
        d12 = np.column_stack([self.x12, zero, zero])
        d43 = np.column_stack([self.x43, zero, zero])
        coord_id = self.model.aero_coord()
        if coord_id is None:
            # yes, this really does list + array addition
            p2 = p1 + d12
            p3 = p4 + d43
        else:
            coord = self.model.coord
            p2 = p1 + coord.transform_local_xyz_to_global(d12, coord_id)
            p3 = p4 + coord.transform_local_xyz_to_global(d43, coord_id)
        return p1, p2, p3, p4

    @property
    def min_max_eid(self) -> tuple[int, int]:
        """
        Gets the min and max element ids of the CAERO card

        Returns
        -------
        min_max_eid : (2, ) list
            The [min_eid, max_eid]

        """
        nchord, nspan = self.shape
        mins = self.element_id
        maxs = self.element_id + nchord * nspan
        return mins.min(), maxs.max()

    @property
    def shape(self) -> tuple[np.ndarray, np.ndarray]:
        nchord = self.nchord.copy()
        nspan = self.nspan.copy()

        ichord = (nchord == 0)
        ispan = (nspan == 0)

        aefact = self.model.aefact
        if np.any(ichord):
            jchord = np.where(ichord)[0]
            lchord = self.lchord[jchord]
            aefacti = aefact.slice_card_by_aefact_id(lchord)
            nchord[ichord] = aefacti.nfractions - 1  # points -> boxes

        if np.any(ispan):
            jspan = np.where(ispan)[0]
            lspan = self.lspan[jspan]
            aefacti = aefact.slice_card_by_aefact_id(lspan)
            nspan[ispan] = aefacti.nfractions - 1  # points -> boxes
        return nchord, nspan

    def get_npanel_points_elements(self) -> tuple[int, int]:
        nchord, nspan = self.shape

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        assert nelements.min() > 0, nelements
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []
        p1, p2, p3, p4 = self.get_points()

        x, y = self.xy
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        igroup_ = array_str(self.igroup, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nspan_ = array_default_int(self.nspan, default=0, size=size)
        lspan_ = array_default_int(self.lspan, default=0, size=size)
        nchord_ = array_default_int(self.nchord, default=0, size=size)
        lchord_ = array_default_int(self.lchord, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, pid, igroup, p1, x12, p4, x43, cp, nspan, lspan, nchord, lchord in zip(\
                element_id, property_id, igroup_, p1_, self.x12, p4_, self.x43, cp_,
                nspan_, lspan_, nchord_, lchord_):
            list_fields = ['CAERO1', eid, pid, cp, nspan, nchord,
                           lspan, lchord, igroup] + p1 + [x12] + p4 + [x43]
            bdf_file.write(print_card(list_fields))
        return


class CAERO2(VectorizedBaseCard):
    """
    Aerodynamic Body Connection
    Defines aerodynamic slender body and interference elements for
    Doublet-Lattice aerodynamics.

    +--------+-----+-----+----+-----+------+-----+------+------+
    |    1   |  2  |  3  |  4 |  5  |   6  |   7 |   8  |  9   |
    +========+=====+=====+====+=====+======+=====+======+======+
    | CAERO2 | EID | PID | CP | NSB | NINT | LSB | LINT | IGID |
    +--------+-----+-----+----+-----+------+-----+------+------+
    |        | X1  |  Y1 | Z1 | X12 |      |     |      |      |
    +--------+-----+-----+----+-----+------+-----+------+------+
    """
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.igroup = np.array([], dtype='int32')
        self.p1 = np.zeros((0, 3), dtype='float64')
        self.x12 = np.array([], dtype='float64')
        self.cp = np.array([], dtype='int32')

        self.nsb = np.array([], dtype='int32')
        self.nint = np.array([], dtype='int32')
        self.lsb = np.array([], dtype='int32')
        self.lint = np.array([], dtype='int32')

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nsb = integer_or_blank(card, 4, 'nsb', default=0)
        nint = integer_or_blank(card, 5, 'nint', default=0)

        lsb = integer_or_blank(card, 6, 'nsb=%s lsb' % nsb, default=0)
        lint = integer_or_blank(card, 7, 'nint=%s lint' % nint, default=0)
        igroup = integer(card, 8, 'igroup')

        p1 = np.array([
            double_or_blank(card, 9, 'x1', default=0.0),
            double_or_blank(card, 10, 'y1', default=0.0),
            double_or_blank(card, 11, 'z1', default=0.0)])
        x12 = double_or_blank(card, 12, 'x12', default=0.)
        assert len(card) <= 13, f'len(CAERO2 card) = {len(card):d}\ncard={card}'
        #return CAERO2(eid, pid, igroup, p1, x12,
                      #cp=cp, nsb=nsb, nint=nint, lsb=lsb, lint=lint,
                      #comment=comment)

        self.cards.append((eid, pid, igroup, p1, x12,
                           cp, nsb, nint, lsb, lint, comment))
        self.n += 1
        return self.n - 1

    def add(self, eid: int, pid: int, igroup: int,
            p1: list[float], x12: float,
            cp: int=0,
            nsb: int=0, nint: int=0,
            lsb: int=0, lint: int=0, comment: str='') -> int:
        """
        Defines a CAERO2 card, which defines a slender body
        (e.g., fuselage/wingtip tank).

        Parameters
        ----------
        eid : int
            element id
        pid : int, PAERO2
            int : PAERO2 ID
        igroup : int
            Group number
        p1 : (1, 3) ndarray float
            xyz location of point 1 (forward position)
        x12 : float
            length of the CAERO2
        cp : int; default=0
            int : coordinate system
        nsb : int; default=0
            Number of slender body elements
        lsb : int; default=0
            AEFACT id for defining the location of the slender body elements
        nint : int; default=0
            Number of interference elements
        lint : int; default=0
            AEFACT id for defining the location of interference elements
        comment : str; default=''
            a comment for the card

        """
        nsb = nsb if nsb is not None else 0
        nint = nint if nint is not None else 0
        lsb = lsb if lsb is not None else 0
        lint = lint if lint is not None else 0
        self.cards.append((eid, pid, igroup, p1, x12,
                           cp, nsb, nint, lsb, lint, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        igroup = np.zeros(ncards, dtype='int32')
        p1 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nsb = np.zeros(ncards, dtype='int32')
        nint = np.zeros(ncards, dtype='int32')
        lsb = np.zeros(ncards, dtype='int32')
        lint = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (eid, pid, igroupi, p1i, x12i,
             cpi, nsbi, ninti, lsbi, linti, comment) = card
            element_id[icard] = eid
            property_id[icard] = pid
            igroup[icard] = igroupi
            p1[icard, :] = p1i
            x12[icard] = x12i

            cp[icard] = cpi
            nsb[icard] = nsbi
            nint[icard] = ninti
            lsb[icard] = lsbi
            lint[icard] = linti
        self._save(element_id, property_id, igroup, p1, x12,
                   cp, nsb, nint, lsb, lint)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, igroup, p1, x12,
              cp, nsb, nint, lsb, lint) -> None:
        if len(self.element_id):
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            igroup = np.hstack([self.igroup, igroup])
            p1 = np.vstack([self.p1, p1])
            x12 = np.hstack([self.x12, x12])
            cp = np.hstack([self.cp, cp])
            nsb = np.hstack([self.nsb, nsb])
            lsb = np.hstack([self.lsb, lsb])
            nint = np.hstack([self.nint, nint])
            lint = np.hstack([self.lint, lint])

        self.element_id = element_id
        self.property_id = property_id
        self.igroup = igroup
        self.p1 = p1
        self.x12 = x12
        self.cp = cp
        self.nsb = nsb
        self.lsb = lsb
        self.nint = nint
        self.lint = lint
        self.n = len(element_id)

    def __apply_slice__(self, elem: CAERO2, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.igroup = self.igroup[i]
        elem.p1 = self.p1[i, :]
        elem.x12 = self.x12[i]
        elem.cp = self.cp[i]
        elem.nsb = self.nsb[i]
        elem.lsb = self.lsb[i]
        elem.nint = self.nint[i]
        elem.lint = self.lint[i]

    def validate(self) -> None:
        #print('nsb=%s lsb=%s' % (self.nsb, self.lsb))
        #print('nint=%s lint=%s' % (self.nint, self.lint))
        #assert isinstance(self.lsb, integer_types), self.lsb
        #assert isinstance(self.lint, integer_types), self.lint
        #assert len(self.p1) == 3, 'CAERO2: p1=%s' % self.p1
        nsb_lsb = (self.nsb == 0) & (self.lsb == 0)
        nint_lint = (self.nint == 0) & (self.lint == 0)

        if np.any(nsb_lsb):
            msg = 'CAERO2: eid=%s nsb=%s lsb=%s; nsb or lsb must be > 0' % (self.element_id[nsb_lsb], self.nsb[nsb_lsb], self.lsb[nsb_lsb])
            raise ValueError(msg)
        if np.any(nint_lint):
            msg = 'CAERO2: eid=%s nint=%s lint=%s; nint or lint must be > 0' % (self.element_id[nint_lint], self.nint[nint_lint], self.lint[nint_lint])
            raise ValueError(msg)
        #assert len(self.p1) == 3, 'CAERO2: p1=%s' % self.p1
        igroup = (self.igroup < 0)
        if np.any(igroup):
            msg = 'CAERO2: eid=%s nint=%s lint=%s; nint or lint must be > 0' % (self.element_id[igroup], self.igroup[igroup])
            raise ValueError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def get_points(self):
        """creates a 1D representation of the CAERO2"""
        icp = (self.cp != 0)
        if np.any(icp):
            coord = self.model.coord
            ucp = np.unique(self.cp[icp])
            if len(ucp) == 1:
                p1 = coord.transform_local_xyz_to_global(self.p1[icp, :], ucp[0])
            else:
                p1 = coord.transform_local_xyz_to_global(self.p1[icp, :], self.cp[icp])
        else:
            p1 = self.p1
        #p2 = p1.copy()
        #p2[:, 0] = p2[:, 0] + self.x12

        coord_id = self.model.aero_coord()

        d12 = np.zeros(p1.shape)
        d12[:, 0] = self.x12

        p2 = p1 + coord.transform_local_xyz_to_global(d12, coord_id)
        #p2 = p1 + self.ascid_ref.transform_vector_to_global(np.array([self.x12, 0., 0.]))

        #print("x12 = %s" % self.x12)
        #print("pcaero[%s] = %s" % (self.eid, [p1,p2]))
        return (p1, p2)

    def get_points_elements_3d(self):
        """
        Gets the points/elements in 3d space as CQUAD4s
        The idea is that this is used by the GUI to display CAERO panels.

        TODO: doesn't support the aero coordinate system

        """
        paero2s = self.model.paero2.slice_card_by_property_id(self.property_id)  # type: PAERO2
        aefact = self.model.aefact
        p1, p2 = self.get_points()
        L = p2 - p1

        for i, nsb, lsb in zip(count(), self.nsb, self.lsb):
            if nsb == 0:
                aefacti = aefact.slice_card_by_aefact_id(lsb)
                station = aefacti.fractions
                nsb = len(station) - 1
                #print('xstation = ', xstation)
            else:
                station = np.linspace(0., nsb, num=nsb+1) # *dx?
            assert nsb > 0, 'nsb=%s' % nsb

            #print('paero2 - pid=%s lrsb=%s lrib=%s' % (paero2.pid, paero2.lrsb, paero2.lrib))
            #paero2 = paero2s[ipid]
            if paero2s.lrsb[i] == 0:
                radii_slender = np.ones(nsb + 1) * paero2s.width[i]
            else:
                lrsb = paero2s.lrsb[i]
                lrsb_ref = aefact.slice_card_by_aefact_id(lrsb)
                radii_slender = lrsb_ref.fractions

            # TODO: not supported
            if paero2s.lrib[i] == 0:
                unused_radii_interference = np.ones(nsb + 1) * paero2s.width[i]
            else:
                #print('lrib = ', paero2.lrib)
                lrib = paero2s.lrib[i]
                lrib_ref = aefact.slice_card_by_aefact_id(lrib)
                unused_radii_interference = lrib_ref.fractions
            radii = radii_slender

            # TODO: not supported
            #theta_interference1 = paero2.theta1
            #theta_interference2 = paero2.theta2

            if nsb != 0:
                #print('L=%s nx=%s' % (L, nx))
                dxyz = L[i] / nsb
                #print('dxyz\n%s' % (dxyz))
                dx, dy, dz = dxyz
                xstation = station * dx
                ystation = station * dy
                zstation = station * dz
            else:
                dxi = xstation.max() - xstation.min()

                #print('L=%s nx=%s dxi=%s' % (L, nx, dxi))
                xratio = xstation / dxi
                #print('xstation/dxi=%s' % xratio)
                dxyz = np.zeros((nsb+1, 3))
                for i, xr in enumerate(xratio):
                    dxyz[i, :] = xr * L
                ystation = dxyz[:, 1]
                zstation = dxyz[:, 2]

            # I think this just lets you know what directions it can pivot in
            # and therefore doesn't affect visualization
            #assert paero2.orient == 'ZY', paero2.orient
            aspect_ratio = paero2s.aspect_ratio[i]

            assert len(radii) == (nsb + 1), 'len(radii)=%s nx=nsb=%s' % (len(radii), nsb)
            if len(xstation) != (nsb + 1):
                msg = 'len(xstation)=%s nx=nsb=%s\nxstation=%s\n%s' % (
                    len(xstation), nsb, xstation, str(self))
                raise RuntimeError(msg)

            xyz, elems = create_axisymmetric_body(
                xstation, ystation, zstation, radii, aspect_ratio,
                p1[i, :])
            assert xyz is not None, str(self)
        return xyz, elems

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        igroup_ = array_str(self.igroup, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nsb_ = array_default_int(self.nsb, default=0, size=size)
        lsb_ = array_default_int(self.lsb, default=0, size=size)
        nint_ = array_default_int(self.nint, default=0, size=size)
        lint_ = array_default_int(self.lint, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        p1_ = self.p1.tolist()
        for eid, pid, igroup, p1, x12, cp, nsb, lsb, nint, lint in zip(\
            element_id, property_id, igroup_, p1_, self.x12, cp_,
            nsb_, lsb_, nint_, lint_):
            list_fields = (['CAERO2', eid, pid, cp, nsb, nint,
                            lsb, lint, igroup, ] + list(p1) +
                           [x12])
            bdf_file.write(print_card(list_fields))
        return


class CAERO3(VectorizedBaseCard):
    """
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO4 | EID | PID | CP | LISTW | LISTC1 | LISTC2 |        |      |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int
        int : PAERO1 ID
    igroup : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int
        int : coordinate system
    list_w : int
        Identification number of an AEFACT entry that lists (x,y) pairs for
        structural interpolation of the wing. (Integer > 0)
    LISTC1, LISTC2 : int
        Identification number of AEFACT entries that list (x,y) pairs for
        control surfaces, if they exist. (Integer >= 0)
    comment : str; default=''
         a comment for the card

    """
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')
        self.property_id = np.array([], dtype='int32')
        self.p1 = np.zeros((0, 3), dtype='float64')
        self.p4 = np.zeros((0, 3), dtype='float64')
        self.x12 = np.array([], dtype='float64')
        self.x43 = np.array([], dtype='float64')
        self.cp = np.array([], dtype='int32')

        self.list_w = np.array([], dtype='int32')
        self.list_c1 = np.array([], dtype='int32')
        self.list_c2 = np.array([], dtype='int32')

    def add(self, eid: int, pid: int,
            p1: np.ndarray, x12: float,
            p4: np.ndarray, x43: float,
            cp: int=0, list_w: int=0, list_c1=None, list_c2=None, comment='') -> int:
        """Creates a CAERO3 card"""
        list_w = list_c1 if list_w else 0
        list_c1 = list_c1 if list_c1 else 0
        list_c2 = list_c2 if list_c2 else 0
        assert isinstance(list_w, integer_types), list_w
        assert isinstance(list_c1, integer_types), list_c1
        assert isinstance(list_c2, integer_types), list_c2
        card = (eid, pid, p1, x12, p4, x43, cp, list_w, list_c1, list_c2, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', default=0)

        list_w = integer_or_blank(card, 4, 'list_w', default=0)
        list_c1 = integer_or_blank(card, 5, 'list_c1', default=0)  ## TODO: should this have a default?
        list_c2 = integer_or_blank(card, 6, 'list_c2', default=0)  ## TODO: should this have a default?

        p1 = np.array([
            double_or_blank(card, 9, 'x1', default=0.0),
            double_or_blank(card, 10, 'y1', default=0.0),
            double_or_blank(card, 11, 'z1', default=0.0)])
        x12 = double_or_blank(card, 12, 'x12', default=0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', default=0.0),
            double_or_blank(card, 14, 'y4', default=0.0),
            double_or_blank(card, 15, 'z4', default=0.0)])
        x43 = double_or_blank(card, 16, 'x43', default=0.)

        assert len(card) <= 17, f'len(CAERO3 card) = {len(card):d}\ncard={card}'
        card = (eid, pid, p1, x12, p4, x43, cp, list_w, list_c1, list_c2, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        list_w = np.zeros(ncards, dtype='int32')
        list_c1 = np.zeros(ncards, dtype='int32')
        list_c2 = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            eid, pid, p1i, x12i, p4i, x43i, cpi, list_wi, list_c1i, list_c2i, comment = card
            element_id[icard] = eid
            property_id[icard] = pid
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            list_w[icard] = list_wi
            list_c1[icard] = list_c1i
            list_c2[icard] = list_c2i
        self._save(element_id, property_id, p1, p4, x12, x43, cp,
                   list_w, list_c1, list_c2)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, p1, p4, x12, x43, cp,
              list_w, list_c1, list_c2):
        if len(self.element_id) != 0:
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            p1 = np.vstack([self.p1, p1])
            p4 = np.vstack([self.p4, p4])
            x12 = np.hstack([self.x12, x12])
            x43 = np.hstack([self.x43, x43])
            cp = np.hstack([self.cp, cp])
            list_w = np.hstack([self.list_w, list_w])
            list_c1 = np.hstack([self.list_c1, list_c1])
            list_c2 = np.hstack([self.list_c2, list_c2])

        self.element_id = element_id
        self.property_id = property_id
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.list_w = list_w
        self.list_c1 = list_c1
        self.list_c2 = list_c2
        self.n = len(element_id)

    def __apply_slice__(self, elem: CAERO3, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.list_w = self.list_w[i]
        elem.list_c1 = self.list_c1[i]
        elem.list_c2 = self.list_c2[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def get_leading_edge_points(self):
        """gets the leading edge points"""
        p1 = self.p1
        p4 = self.p4
        icps = (self.cp != 0)
        if np.any(icps):
            cps = self.cp[icps]
            coord = self.model.coord
            for icp, cp in zip(icps, cps):
                #coord = coord.slice_card_by_coord_id(self.cp)
                p1[icp, :] = coord.transform_local_xyz_to_global(self.p1[icp, :], cp)
                p4[icp, :] = coord.transform_local_xyz_to_global(self.p4[icp, :], cp)
        return p1, p4

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        p1, p4 = self.get_leading_edge_points()
        zero = np.zeros(self.x12.size, dtype=self.x12.dtype)
        d12 = np.column_stack([self.x12, zero, zero])
        d43 = np.column_stack([self.x43, zero, zero])
        coord_id = self.model.aero_coord()
        if coord_id is None:
            # yes, this really does list + array addition
            p2 = p1 + d12
            p3 = p4 + d43
        else:
            coord = self.model.coord
            p2 = p1 + coord.transform_local_xyz_to_global(d12, coord_id)
            p3 = p4 + coord.transform_local_xyz_to_global(d43, coord_id)
        return p1, p2, p3, p4

    @property
    def xy(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        nchord, nspan = self.shape
        x = []
        y = []
        xi = np.linspace(0., 1., nchord + 1)
        for eid, nspani in zip(self.element_id, nspan):
            yi = np.linspace(0., 1., nspani + 1)
            if nspani < 1:
                msg = 'CAERO3 eid=%s nspan=%s' % (eid, nspani)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    @property
    def shape(self) -> tuple[int, np.ndarray]:
        nchord = 2
        paero3 = self.model.paero3.slice_card_by_property_id(self.property_id)  # type: PAERO3
        nspan = paero3.nbox
        return nchord, nspan

    def get_npanel_points_elements(self) -> tuple[int, int]:
        nchord, nspan = self.shape

        #aefact = self.model.aefact
        #if np.any(ichord):
            #jchord = np.where(ichord)[0]
            #lchord = self.lchord[jchord]
            #aefacti = aefact.slice_card_by_aefact_id(lchord)
            #nchord[ichord] = aefacti.nfractions - 1  # points -> boxes

            #for jchordi, lchordi in zip(jchord, lchord):

        #if np.any(ispan):
            #jspan = np.where(ispan)[0]
            #lspan = self.lspan[jspan]
            #aefacti = aefact.slice_card_by_aefact_id(lspan)
            #nspan[ispan] = aefacti.nfractions - 1  # points -> boxes

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        assert nelements.min() > 0, nelements
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []
        p1, p2, p3, p4 = self.get_points()

        x, y = self.xy
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        list_w_ = array_default_int(self.list_w, default=0, size=size)
        list_c1_ = array_default_int(self.list_c1, default=0, size=size)
        list_c2_ = array_default_int(self.list_c2, default=0, size=size)

        #nspan_ = array_default_int(self.nspan, default=0, size=size)
        #lspan_ = array_default_int(self.lspan, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, pid, p1, x12, p4, x43, cp, list_w, list_c1, list_c2 in zip(
                element_id, property_id, p1_, self.x12, p4_, self.x43, cp_,
                list_w_, list_c1_, list_c2_):
            list_fields = ['CAERO3', eid, pid, cp, list_w, list_c1, list_c2,
                           None, None] + p1 + [x12] + p4 + [x43]
            #list_fields = ['CAERO3', eid, pid, cp, nspan, lspan,
                           #'', '', ''] + p1 + [x12] + p4 + [x43]
            bdf_file.write(print_card(list_fields))
        return


class CAERO4(VectorizedBaseCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO4 | EID | PID | CP | NSPAN |  LSPAN |        |        |      |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int
        int : PAERO1 ID
    igroup : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int
        int : coordinate system
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
    comment : str; default=''
         a comment for the card

    """
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int,
            p1: np.ndarray, x12: float,
            p4: np.ndarray, x43: float,
            cp: int=0, nspan: int=0, lspan: int=0, comment: str='') -> int:
        """
        Defines a CAERO4 card, which defines a strip theory surface.

        Parameters
        ----------
        eid : int
            element id
        pid : int
            int : PAERO4 ID
        p1 : (1, 3) ndarray float
            xyz location of point 1 (leading edge; inboard)
        p4 : (1, 3) ndarray float
            xyz location of point 4 (leading edge; outboard)
        x12 : float
            distance along the flow direction from node 1 to node 2
            (typically x, root chord)
        x43 : float
            distance along the flow direction from node 4 to node 3
            (typically x, tip chord)
        cp : int; default=0
            int : coordinate system
        nspan : int; default=0
            int > 0 : N spanwise boxes distributed evenly
            int = 0 : use lchord
        lspan : int; default=0
            int > 0 : AEFACT reference for non-uniform nspan
            int = 0 : use nspan
        comment : str; default=''
             a comment for the card

        """
        card = (eid, pid, p1, x12, p4, x43, cp, nspan, lspan, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        lspan = integer_or_blank(card, 5, 'lspan', default=0)

        p1 = np.array([
            double_or_blank(card, 9, 'x1', default=0.0),
            double_or_blank(card, 10, 'y1', default=0.0),
            double_or_blank(card, 11, 'z1', default=0.0)])
        x12 = double_or_blank(card, 12, 'x12', default=0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', default=0.0),
            double_or_blank(card, 14, 'y4', default=0.0),
            double_or_blank(card, 15, 'z4', default=0.0)])
        x43 = double_or_blank(card, 16, 'x43', default=0.)

        assert len(card) <= 17, f'len(CAERO4 card) = {len(card):d}\ncard={card}'
        #return CAERO1(eid, pid, igroup, p1, x12, p4, x43,
                      #cp=cp, nspan=nspan, lspan=lspan, nchord=nchord, lchord=lchord,
                      #comment=comment)
        card = (eid, pid, p1, x12, p4, x43, cp, nspan, lspan, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nspan = np.zeros(ncards, dtype='int32')
        lspan = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            eid, pid, p1i, x12i, p4i, x43i, cpi, nspani, lspani, comment = card
            element_id[icard] = eid
            property_id[icard] = pid
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            nspan[icard] = nspani
            lspan[icard] = lspani
        self._save(element_id, property_id, p1, p4, x12, x43, cp,
                   nspan, lspan)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, p1, p4, x12, x43, cp,
              nspan, lspan):
        if len(self.element_id):
            element_id = np.hstack([self.element_id, element_id])
            property_id = np.hstack([self.property_id, property_id])
            p1 = np.vstack([self.p1, p1])
            p4 = np.vstack([self.p4, p4])
            x12 = np.hstack([self.x12, x12])
            x43 = np.hstack([self.x43, x43])
            cp = np.hstack([self.cp, cp])
            nspan = np.hstack([self.nspan, nspan])
            lspan = np.hstack([self.lspan, lspan])

        self.element_id = element_id
        self.property_id = property_id
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.n = len(element_id)

    def __apply_slice__(self, elem: CAERO4, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.nspan = self.nspan[i]
        elem.lspan = self.lspan[i]

    def validate(self):
        msg = ''
        is_failed = False
        #if not isinstance(self.p1, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True
        #if not isinstance(self.p4, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True

        element_id = self.element_id

        ibad = (self.x12 <= 0.)
        if np.any(ibad):
            msg += 'X12 and must be greater than or equal to 0\n'
            msg += f'  element_id={element_id[ibad]}\n   x12={self.x12[ibad]}'
            is_failed = True
        #if self.x43 <= 0.:
            #msg += 'X43=%s and must be greater than or equal to 0\n' % (self.x43)
            #is_failed = True

        ibad = (self.nspan == 0) & (self.lspan == 0)
        if np.any(ibad):
            msg += 'NSPAN or LSPAN must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nspan != 0) & (self.lspan != 0)
        if np.any(ibad):
            msg += 'Either NSPAN or LSPAN must 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        if is_failed:
            msg += str(self)
            #msg += CAERO4_MSG
            raise RuntimeError(msg)

        assert self.p1.shape[1] == 3, 'p1=%s' % self.p1.shape
        assert self.p4.shape[1] == 3, 'p4=%s' % self.p4.shape

        # calculating area; assuming coordinate transformations don't matter
        if 1:
            p1 = self.p1
            p4 = self.p4
            d12 = np.zeros(p1.shape)
            d12[:, 0] = self.x12
            d43 = np.zeros(p1.shape)
            d43[:, 0] = self.x43
            p2 = p1 + d12
            p3 = p4 + d43

            a = p3 - p1
            b = p4 - p2
            area = np.linalg.norm(np.cross(a, b), axis=1)
            assert len(area) == p1.shape[0]
            ibad = (area < 0.0)
            if np.any(ibad):
                msg += 'Either NCHORD or LCHORD must 0\n'
                msg += f'eid={self.element_id[ibad,:]} p1={p1[ibad,:]} p2={p2[ibad,:]} p3={p3[ibad,:]} p4={p4[ibad,:]} area={area[ibad]}\n'
                raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    @property
    def shape(self) -> tuple[np.ndarray, np.ndarray]:
        nchord = 1
        nspan = self.nspan.copy()

        #ichord = (nchord == 0)
        ispan = (nspan == 0)

        aefact = self.model.aefact
        if np.any(ispan):
            jspan = np.where(ispan)[0]
            lspan = self.lspan[jspan]
            aefacti = aefact.slice_card_by_aefact_id(lspan)
            nspan[ispan] = aefacti.nfractions - 1  # points -> boxes
        return nchord, nspan

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        x = []
        y = []
        aefact = self.model.aefact
        nchord = 1
        xi = np.linspace(0., 1., nchord + 1)  # nchord=1
        for eid, nspan, lspan in zip(self.element_id, self.nspan, self.lspan):
            if nspan == 0:
                aefacti = aefact.slice_card_by_aefact_id(lspan)
                nspan = aefacti.nfractions - 1  # points -> boxes
                yi = aefacti.fractions
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nspan < 1:
                msg = 'CAERO4 eid=%s nspan=%s lspan=%s' % (eid, nspan, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    def get_leading_edge_points(self):
        """gets the leading edge points"""
        p1 = self.p1
        p4 = self.p4
        icps = (self.cp != 0)
        if np.any(icps):
            #coord = self.model.coord.slice_card_by_coord_id(self.cp)
            cps = self.cp[icps]
            for icp, cp in zip(icps, cps):
                p1[icp, :] = self.model.coord.transform_local_xyz_to_global(self.p1[icp, :], cp)
                p4[icp, :] = self.model.coord.transform_local_xyz_to_global(self.p4[icp, :], cp)
        return p1, p4

    def get_points(self):
        """
        Get the 4 corner points for the CAERO card

        Returns
        -------
        p1234 : (4, 3) list
             List of 4 corner points in the global frame

        """
        p1, p4 = self.get_leading_edge_points()
        zero = np.zeros(self.x12.size, dtype=self.x12.dtype)
        d12 = np.column_stack([self.x12, zero, zero])
        d43 = np.column_stack([self.x43, zero, zero])
        coord_id = self.model.aero_coord()
        if coord_id is None:
            # yes, this really does list + array addition
            p2 = p1 + d12
            p3 = p4 + d43
        else:
            coord = self.model.coord
            p2 = p1 + coord.transform_local_xyz_to_global(d12, coord_id)
            p3 = p4 + coord.transform_local_xyz_to_global(d43, coord_id)
        return p1, p2, p3, p4

    def get_npanel_points_elements(self) -> tuple[int, int]:
        nchord, nspan = self.shape

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        assert nelements.min() > 0, nelements
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []
        p1, p2, p3, p4 = self.get_points()

        x, y = self.xy
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nspan_ = array_default_int(self.nspan, default=0, size=size)
        lspan_ = array_default_int(self.lspan, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, pid, p1, x12, p4, x43, cp, nspan, lspan in zip(
                element_id, property_id, p1_, self.x12, p4_, self.x43, cp_,
                nspan_, lspan_):
            list_fields = ['CAERO4', eid, pid, cp, nspan, lspan,
                           '', '', ''] + p1 + [x12] + p4 + [x43]
            bdf_file.write(print_card(list_fields))
        return


class CAERO5(VectorizedBaseCard):
    """
    Defines an aerodynamic macro element (panel) in terms of two leading edge
    locations and side chords. This is used for Doublet-Lattice theory for
    subsonic aerodynamics and the ZONA51 theory for supersonic aerodynamics.

    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |   1    |  2  |  3  | 4  |   5   |   6    |    7   |   8    |   9  |
    +========+=====+=====+====+=======+========+========+========+======+
    | CAERO5 | EID | PID | CP | NSPAN | LSPAN  |  NTHRY | NTHICK |      |
    +--------+-----+-----+----+-------+--------+--------+--------+------+
    |        |  X1 | Y1  | Z1 |  X12  |   X4   |   Y4   |   Z4   | X43  |
    +--------+-----+-----+----+-------+--------+--------+--------+------+

    ::

      1
      | \
      |   \
      |     \
      |      4
      |      |
      |      |
      2------3

    Attributes
    ----------
    eid : int
        element id
    pid : int
        int : PAERO1 ID
    igroup : int
        Group number
    p1 : (1, 3) ndarray float
        xyz location of point 1 (leading edge; inboard)
    p4 : (1, 3) ndarray float
        xyz location of point 4 (leading edge; outboard)
    x12 : float
        distance along the flow direction from node 1 to node 2; (typically x, root chord)
    x43 : float
        distance along the flow direction from node 4 to node 3; (typically x, tip chord)
    cp : int
        int : coordinate system
    nspan : int
        int > 0 : N spanwise boxes distributed evenly
        int = 0 : use lchord
    lspan : int
        int > 0 : AEFACT reference for non-uniform nspan
        int = 0 : use nspan
    comment : str; default=''
         a comment for the card

    """
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')

    def add(self, eid: int, pid: int,
            p1: list[float], x12: float,
            p4: list[float], x43: float,
            cp: int=0,
            nspan: int=0, lspan: int=0,
            ntheory: int=0,
            nthick: int=0, comment: str='') -> int:
        """Creates a CAERO5 card"""
        card = (eid, pid, p1, x12, p4, x43, cp, nspan, lspan,
                ntheory, nthick, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        pid = integer(card, 2, 'pid')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        lspan = integer_or_blank(card, 5, 'lspan', default=0)
        ntheory = integer_or_blank(card, 6, 'ntheory', default=0)
        nthick = integer_or_blank(card, 7, 'nthick', default=0)

        p1 = np.array([
            double_or_blank(card, 9, 'x1', default=0.0),
            double_or_blank(card, 10, 'y1', default=0.0),
            double_or_blank(card, 11, 'z1', default=0.0)])
        x12 = double_or_blank(card, 12, 'x12', default=0.)

        p4 = np.array([
            double_or_blank(card, 13, 'x4', 0.0),
            double_or_blank(card, 14, 'y4', 0.0),
            double_or_blank(card, 15, 'z4', 0.0)])
        x43 = double_or_blank(card, 16, 'x43', 0.)

        assert len(card) <= 17, f'len(CAERO5 card) = {len(card):d}\ncard={card}'
        card = (eid, pid, p1, x12, p4, x43, cp, nspan, lspan, ntheory, nthick, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        property_id = np.zeros(ncards, dtype='int32')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nspan = np.zeros(ncards, dtype='int32')
        lspan = np.zeros(ncards, dtype='int32')

        ntheory = np.zeros(ncards, dtype='int32')
        nthick = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (eid, pid, p1i, x12i, p4i, x43i, cpi, nspani, lspani, ntheoryi, nthicki, commenti) = card
            element_id[icard] = eid
            property_id[icard] = pid
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            nspan[icard] = nspani
            lspan[icard] = lspani
            ntheory[icard] = ntheoryi
            nthick[icard] = nthicki
        self._save(element_id, property_id, p1, p4, x12, x43, cp,
                   nspan, lspan, ntheory, nthick)
        self.sort()
        self.cards = []

    def _save(self, element_id, property_id, p1, p4, x12, x43, cp,
              nspan, lspan, ntheory, nthick):
        assert len(self.element_id) == 0, self.element_id
        self.element_id = element_id
        self.property_id = property_id
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.ntheory = ntheory
        self.nthick = nthick
        self.n = len(element_id)

    def __apply_slice__(self, elem: CAERO5, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.property_id[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.nspan = self.nspan[i]
        elem.lspan = self.lspan[i]
        elem.ntheory = self.ntheory[i]
        elem.nthick = self.nthick[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        x = []
        y = []
        aefact = self.model.aefact
        nchord = 1
        xi = np.linspace(0., 1., nchord + 1)

        for eid, nspan, lspan in zip(self.element_id, self.nspan, self.lspan):
            if nspan == 0:
                aefacti = aefact.slice_card_by_aefact_id(lspan)
                nspan = aefacti.nfractions - 1  # points -> boxes
                yi = aefacti.fractions
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nchord < 1 or nspan < 1:
                msg = 'CAERO5 eid=%s nchord=%s nspan=%s lspan=%s' % (
                    eid, nchord, nspan, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    def get_npanel_points_elements(self) -> tuple[int, int]:
        #nchord = self.nchord.copy()
        nspan = self.nspan.copy()
        nchord = np.ones(nspan.size, dtype=nspan.dtype)

        #ichord = (nchord == 0)
        ispan = (nspan == 0)

        aefact = self.model.aefact
        #if np.any(ichord):
            #jchord = np.where(ichord)[0]
            #lchord = self.lchord[jchord]
            #aefacti = aefact.slice_card_by_aefact_id(lchord)
            #nchord[ichord] = aefacti.nfractions - 1  # points -> boxes

        if np.any(ispan):
            jspan = np.where(ispan)[0]
            lspan = self.lspan[jspan]
            aefacti = aefact.slice_card_by_aefact_id(lspan)
            nspan[ispan] = aefacti.nfractions - 1  # points -> boxes

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        assert nelements.min() > 0, nelements
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []

        icp = (self.cp != 0)
        if np.any(icp):
            coord = self.model.coord
            raise NotImplementedError('CAERO1: cp')
        else:
            p1 = self.p1
            p4 = self.p4
            p2 = p1.copy()
            p3 = p4.copy()
        p2[:, 0] = p2[:, 0] + self.x12
        p4[:, 0] = p4[:, 0] + self.x43

        x, y = self.xy
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        property_id = array_str(self.property_id, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nspan_ = array_default_int(self.nspan, default=0, size=size)
        lspan_ = array_default_int(self.lspan, default=0, size=size)
        nthick_ = array_default_int(self.ntheory, default=0, size=size)
        ntheory_ = array_default_int(self.nthick, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, pid, p1, x12, p4, x43, cp, nspan, lspan, ntheory, nthick in zip(\
                element_id, property_id, p1_, self.x12, p4_, self.x43, cp_,
                nspan_, lspan_, ntheory_, nthick_):
            list_fields = ['CAERO5', eid, pid, cp, nspan, lspan,
                           ntheory, nthick, ''] + p1 + [x12] + p4 + [x43]
            bdf_file.write(print_card(list_fields))
        return


class CAERO7(VectorizedBaseCard):
    _id_name = 'element_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.element_id = np.array([], dtype='int32')

    def add(self, eid: int, label: str,
            p1: np.ndarray, x12: float,
            p4: np.ndarray, x43: float,
            cp: int=0,
            nspan: int=0, lspan: int=0,
            nchord: int=0,
            p_airfoil: int=0, ztaic: int=0, comment: str='') -> int:
        p_airfoil = p_airfoil if p_airfoil is not None else 0
        ztaic = ztaic if ztaic is not None else 0
        card = (eid, label, p1, x12, p4, x43, cp,
                nspan, lspan, nchord,
                p_airfoil, ztaic, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CAERO7 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        eid = integer(card, 1, 'eid')
        name = string(card, 2, 'name')
        cp = integer_or_blank(card, 3, 'cp', default=0)
        nspan = integer_or_blank(card, 4, 'nspan', default=0)
        nchord = integer_or_blank(card, 5, 'nchord', default=0)
        lspan = integer_or_blank(card, 6, 'aefact_lchord', default=0)
        ztaic = integer_or_blank(card, 7, 'ztaic', default=0)
        p_airfoil = integer_or_blank(card, 8, 'aefact', default=0)
        #assert cp == 0
        #igroup = integer(card, 8, 'igid')

        x1 = double_or_blank(card, 9, 'x1', default=0.0)
        y1 = double_or_blank(card, 10, 'y1', default=0.0)
        z1 = double_or_blank(card, 11, 'z1', default=0.0)
        p1 = np.array([x1, y1, z1])
        x12 = double_or_blank(card, 12, 'x12', default=0.)
        unused_lchord_root = integer_or_blank(card, 13, 'lchord_root')
        unused_attach_root = integer_or_blank(card, 14, 'attach_root')
        unused_achord_root = integer_or_blank(card, 15, 'achord_root')

        x4 = double_or_blank(card, 17, 'x4', default=0.0)
        y4 = double_or_blank(card, 18, 'y4', default=0.0)
        z4 = double_or_blank(card, 19, 'z4', default=0.0)
        p4 = np.array([x4, y4, z4])
        x43 = double_or_blank(card, 20, 'x43', default=0.)
        unused_lchord_tip = integer_or_blank(card, 21, 'lchord_tip')
        unused_attach_tip = integer_or_blank(card, 22, 'attach_tip')
        unused_achord_tip = integer_or_blank(card, 23, 'achord_tip')

        assert len(card) <= 23, f'len(CAERO7 card) = {len(card):d}\ncard={card}'
        #return CAERO7(eid, name, p1, x12, p4, x43,
                      #cp=cp, nspan=nspan, nchord=nchord, lspan=lspan,
                      #p_airfoil=p_airfoil, ztaic=ztaic,
                      #comment=comment)
        card = (eid, name, p1, x12, p4, x43,
                cp, nspan, lspan, nchord,
                p_airfoil, ztaic, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        element_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        p1 = np.zeros((ncards, 3), dtype='float64')
        p4 = np.zeros((ncards, 3), dtype='float64')
        x12 = np.zeros(ncards, dtype='float64')
        x43 = np.zeros(ncards, dtype='float64')
        cp = np.zeros(ncards, dtype='int32')

        nspan = np.zeros(ncards, dtype='int32')
        nchord = np.zeros(ncards, dtype='int32')
        lspan = np.zeros(ncards, dtype='int32')
        lchord = np.zeros(ncards, dtype='int32')
        ztaic = np.zeros(ncards, dtype='int32')
        p_airfoil = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (eid, labeli, p1i, x12i, p4i, x43i,
             cpi, nspani, lspani, nchordi,
             p_airfoili, ztaici, comment) = card
            assert len(labeli) <= 8, f'label={labeli!r}'
            element_id[icard] = eid
            label[icard] = labeli
            p1[icard, :] = p1i
            x12[icard] = x12i
            p4[icard, :] = p4i
            x43[icard] = x43i
            cp[icard] = cpi
            nspan[icard] = nspani
            lspan[icard] = lspani
            nchord[icard] = nchordi
            #lchord[icard] = lchordi
            ztaic[icard] = ztaici
            p_airfoil[icard] = p_airfoili
        self._save(element_id, label, p1, p4, x12, x43, cp,
                   nspan, lspan, nchord, lchord, ztaic, p_airfoil)
        self.sort()
        self.cards = []

    def _save(self, element_id, label, p1, p4, x12, x43, cp,
              nspan, lspan, nchord, lchord, ztaic, p_airfoil):
        assert len(self.element_id) == 0, self.element_id
        self.element_id = element_id
        self.label = label
        self.p1 = p1
        self.p4 = p4
        self.x12 = x12
        self.x43 = x43
        self.cp = cp
        self.nspan = nspan
        self.lspan = lspan
        self.nchord = nchord
        self.lchord = lchord
        self.ztaic = ztaic
        self.p_airfoil = p_airfoil
        self.n = len(element_id)

    def __apply_slice__(self, elem: CAERO7, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.element_id = self.element_id[i]
        elem.property_id = self.label[i]
        elem.p1 = self.p1[i, :]
        elem.p4 = self.p4[i, :]
        elem.x12 = self.x12[i]
        elem.x43 = self.x43[i]
        elem.cp = self.cp[i]
        elem.nspan = self.nspan[i]
        elem.lspan = self.lspan[i]
        elem.nchord = self.nchord[i]
        elem.lchord = self.lchord[i]
        elem.ztaic = self.ztaic[i]
        elem.p_airfoil = self.p_airfoil[i]

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        #xy = []
        x = []
        y = []
        for eid, nchord, lchord, nspan, lspan in zip(self.element_id,
                                                     self.nchord, self.lchord,
                                                     self.nspan, self.lspan):
            if nchord == 0:
                xi = self.lchord_ref.fractions
                nchord = len(xi) - 1
            else:
                xi = np.linspace(0., 1., nchord + 1)

            if nspan == 0:
                yi = self.lspan_ref.fractions
                nspan = len(yi) - 1
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nchord < 1 or nspan < 1:
                msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                    eid, nchord, nspan, lchord, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    @property
    def shape(self) -> tuple[np.ndarray, np.ndarray]:
        nchord = self.nchord.copy()
        nspan = self.nspan.copy()

        ichord = (nchord == 0)
        ispan = (nspan == 0)

        aefact = self.model.aefact
        if np.any(ichord):
            jchord = np.where(ichord)[0]
            lchord = self.lchord[jchord]
            aefacti = aefact.slice_card_by_aefact_id(lchord)
            nchord[ichord] = aefacti.nfractions - 1  # points -> boxes

        if np.any(ispan):
            jspan = np.where(ispan)[0]
            lspan = self.lspan[jspan]
            aefacti = aefact.slice_card_by_aefact_id(lspan)
            nspan[ispan] = aefacti.nfractions - 1  # points -> boxes
        return nchord, nspan

    def get_npanel_points_elements(self) -> tuple[int, int]:
        nchord, nspan = self.shape

        # nchord, nspan are number of boxes in the different directions
        npoints = (nchord + 1) * (nspan + 1)
        nelements = nchord * nspan
        return npoints, nelements

    def panel_points_elements(self) -> tuple[np.ndarray, np.ndarray]:
        points = []
        elements = []

        icp = (self.cp != 0)
        if np.any(icp):
            coord = self.model.coord
            raise NotImplementedError('CAERO7: cp')
        else:
            p1 = self.p1
            p4 = self.p4
            p2 = p1.copy()
            p3 = p4.copy()
        p2[:, 0] = p2[:, 0] + self.x12
        p4[:, 0] = p4[:, 0] + self.x43

        x, y = self.xy
        ipoint = 0
        for p1i, p2i, p3i, p4i, xi, yi in zip(p1, p2, p3, p4, x, y):
            pointsi, elementsi = points_elements_from_quad_points(p1i, p4i, p3i, p2i, yi, xi, dtype='int32')
            points.append(pointsi)
            elements.append(elementsi + ipoint)
            ipoint += len(pointsi)
        return points, elements

    def validate(self):
        msg = ''
        is_failed = False
        #if not isinstance(self.p1, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True
        #if not isinstance(self.p4, np.ndarray):
            #msg += 'p1=%s and must be a numpy array\n' % (self.p1)
            #is_failed = True

        #element_id = self.element_id

        ibad = (self.x12 <= 0.)
        if np.any(ibad):
            msg += 'X12 and must be greater than or equal to 0\n'
            msg += f'  element_id={element_id[ibad]}\n   x12={self.x12[ibad]}'
            is_failed = True
        #if self.x43 <= 0.:
            #msg += 'X43=%s and must be greater than or equal to 0\n' % (self.x43)
            #is_failed = True

        ibad = (self.nspan == 0) & (self.lspan == 0)
        if np.any(ibad):
            msg += 'NSPAN or LSPAN must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nspan != 0) & (self.lspan != 0)
        if np.any(ibad):
            msg += 'Either NSPAN or LSPAN must 0\n'
            msg += f'  element_id={element_id}\n  nspan={self.nspan[ibad]}\n  lspan={self.lspan}\n'
            is_failed = True

        ibad = (self.nchord == 0) & (self.lchord == 0)
        if np.any(ibad):
            msg += 'NCHORD or LCHORD must be greater than 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True

        ibad = (self.nchord != 0) & (self.lchord != 0)
        if np.any(ibad):
            msg += 'Either NCHORD or LCHORD must 0\n'
            msg += f'  element_id={element_id}\n  nchord={self.nchord[ibad]}\n  lchord={self.lchord}\n'
            is_failed = True
        if is_failed:
            msg += str(self)
            #msg += CAERO7_MSG
            raise ValueError(msg)

        #assert self.p1.shape[1] == 3, 'p1=%s' % self.p1.shape
        #assert self.p4.shape[1] == 3, 'p4=%s' % self.p4.shape

        ## calculating area; assuming coordinate transformations don't matter
        #if 1:
            #p1 = self.p1
            #p4 = self.p4
            #d12 = np.zeros(p1.shape)
            #d12[:, 0] = self.x12
            #d43 = np.zeros(p1.shape)
            #d43[:, 0] = self.x43
            #p2 = p1 + d12
            #p3 = p4 + d43

            #a = p3 - p1
            #b = p4 - p2
            #area = np.linalg.norm(np.cross(a, b), axis=1)
            #assert len(area) == p1.shape[0]
            #ibad = (area < 0.0)
            #if np.any(ibad):
                #msg += 'Either NCHORD or LCHORD must 0\n'
                #msg += f'eid={self.element_id[ibad,:]} p1={p1[ibad,:]} p2={p2[ibad,:]} p3={p3[ibad,:]} p4={p4[ibad,:]} area={area[ibad]}\n'
                #raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.cp)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def flip_normal_by_element_id(self, element_id=None):
        """flips the CAERO1 normal vector"""
        i = self.index(element_id)
        self.flip_normal_by_index(i)

    def flip_normal_by_index(self, i: np.ndarray) -> None:
        """flips the CAERO1 normal vector"""
        self.p1[i, :], self.p4[i, :] = self.p4[i, :], self.p1[i, :]
        self.x12[i, :], self.x43[i, :] = self.x43[i, :], self.x12[i, :]

    @property
    def xy(self):
        """
        Returns
        -------
        x : (nchord,) ndarray
            The percentage x location in the chord-wise direction of each panel
        y : (nspan,) ndarray
            The percentage y location in the span-wise direction of each panel

        """
        #xy = []
        x = []
        y = []
        for eid, nchord, lchord, nspan, lspan in zip(self.element_id, self.nchord, self.lchord, self.nspan, self.lspan):
            if nchord == 0:
                xi = self.lchord_ref.fractions
                nchord = len(x) - 1
            else:
                xi = np.linspace(0., 1., nchord + 1)

            if nspan == 0:
                yi = self.lspan_ref.fractions
                nspan = len(y) - 1
            else:
                yi = np.linspace(0., 1., nspan + 1)

            if nchord < 1 or nspan < 1:
                msg = 'CAERO1 eid=%s nchord=%s nspan=%s lchord=%s lspan=%s' % (
                    eid, nchord, nspan, lchord, lspan)
                raise RuntimeError(msg)
            x.append(xi)
            y.append(yi)
        return x, y

    @parse_element_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        element_id = array_str(self.element_id, size=size)
        cp_ = array_default_int(self.cp, default=0, size=size)

        nspan_ = array_default_int(self.nspan, default=0, size=size)
        lspan_ = array_default_int(self.lspan, default=0, size=size)
        nchord_ = array_default_int(self.nchord, default=0, size=size)
        ztaic_ = array_default_int(self.ztaic, default=0, size=size)
        p_airfoil_ = array_default_int(self.p_airfoil, default=0, size=size)
        assert self.p1.shape[1] == 3, self.p1.shape
        assert self.p4.shape[1] == 3, self.p4.shape
        p1_ = self.p1.tolist()
        p4_ = self.p4.tolist()
        for eid, label, p1, x12, p4, x43, cp, nspan, lspan, nchord, ztaic, p_airfoil, in zip(\
                element_id, self.label, p1_, self.x12, p4_, self.x43, cp_,
                nspan_, lspan_, nchord_, ztaic_, p_airfoil_):
            list_fields = [
                'CAERO7', eid, label, cp, nspan, nchord, lspan, ztaic, p_airfoil,] + \
                p1 + [x12, None, None, None, None] + \
                p4 + [x43, None, None, None, None]
            bdf_file.write(print_card(list_fields))
        return


class PAERO(VectorizedBaseCard):
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def slice_card_by_property_id(self, property_id: np.ndarray) -> PAERO:
        assert self.n > 0, self.n
        assert len(self.property_id) > 0, self.property_id
        i = self.index(property_id)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def index(self, property_id: np.ndarray) -> np.ndarray:
        assert len(self.property_id) > 0, self.property_id
        property_id = np.atleast_1d(np.asarray(property_id, dtype=self.property_id.dtype))
        iprop = np.searchsorted(self.property_id, property_id)
        return iprop

    @abstractmethod
    def __apply_slice__(self, prop: PAERO, i: np.ndarray) -> None:
        #...
        raise NotImplementedError(f'{self.type}: add __apply_slice__')


class PAERO1(PAERO):
    """
    Defines associated bodies for the panels in the Doublet-Lattice method.

    +--------+-----+----+----+----+----+----+----+
    |    1   |  2  |  3 |  4 |  5 |  6 |  7 |  8 |
    +========+=====+====+====+====+====+====+====+
    | PAERO1 | PID | B1 | B2 | B3 | B4 | B5 | B6 |
    +--------+-----+----+----+----+----+----+----+

    """
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.caero_body_id = np.array([], dtype='int32')

    def add(self, pid: int, caero_body_ids: Optional[list[int]]=None,
            comment: str='') -> int:
        """
        Creates a PAERO1 card, which defines associated bodies for the
        panels in the Doublet-Lattice method.

        Parameters
        ----------
        pid : int
            PAERO1 id
        caero_body_ids : list[int]; default=None
            CAERO2 ids that are within the same IGID group
        comment : str; default=''
            a comment for the card

        """
        #if caero_body_ids is None:
            #caero_body_ids = []
        self.cards.append((pid, caero_body_ids, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        """
        Adds a PAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        caero_body_ids = [interpret_value(field, card) for field in card[2:]]
        caero_body_ids2 = []

        for caero_body_id in caero_body_ids:
            if isinstance(caero_body_id, integer_types) and caero_body_id >= 0:
                caero_body_ids2.append(caero_body_id)
            elif caero_body_id is not None:
                msg = f'invalid caero_body_id value on PAERO1; caero_body_id={caero_body_id!r}'
                raise RuntimeError(msg)
            #else:
                #pass
        #return PAERO1(pid, caero_body_ids, comment=comment)
        self.cards.append((pid, caero_body_ids, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        #caero_body_id = np.zeros(ncards, dtype='int32')
        caero_body_id = []
        ncaero_body_id = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (pid, caero_body_idsi, comment) = card
            if caero_body_idsi is None:
                continue
            ncardsi = len(caero_body_idsi)
            property_id[icard] = pid
            if ncardsi == 0:
                continue

            ncaero_body_id[icard] = len(caero_body_idsi)
            caero_body_id.extend(caero_body_idsi)
        caero_body_id = np.array(caero_body_id, dtype='int32')
        self._save(property_id, caero_body_id, ncaero_body_id)
        self.sort()
        self.cards = []

    def _save(self, property_id, caero_body_id, ncaero_body_id):
        assert len(self.property_id) == 0, self.property_id
        self.property_id = property_id
        self.caero_body_id = caero_body_id
        self.ncaero_body_id = ncaero_body_id
        self.n = len(property_id)

    def __apply_slice__(self, prop: PAERO1, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        icaero_body_id = self.icaero_body_id
        prop.caero_body_id = hslice_by_idim(i, icaero_body_id, self.caero_body_id)
        prop.ncaero_body_id = self.ncaero_body_id[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        #geom_check(
            #missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.element_id, ucaero_ids),
        #)

    @property
    def icaero_body_id(self) -> np.ndarray:
        return make_idim(self.n, self.ncaero_body_id)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        paero_ids = array_str(self.property_id, size=size)
        caero_body_ids = array_str(self.caero_body_id, size=size)

        for pid, (icaero1, icaero2) in zip(paero_ids, self.icaero_body_id):
            if icaero1 != icaero2:
                list_fields = ['PAERO1', pid] + caero_body_ids[icaero1:icaero2].tolist()
            else:
                list_fields = ['PAERO1', pid]
            bdf_file.write(print_card(list_fields))
        return


class PAERO2(PAERO):
    """
    Defines the cross-sectional properties of aerodynamic bodies.

    +--------+------+--------+-------+------+------+------+------+------+
    |    1   |  2   |   3    |    4  |   5  |   6  |   7  |   8  |   9  |
    +========+======+========+=======+======+======+======+======+======+
    | PAERO2 | PID  | ORIENT | WIDTH |  AR  | LRSB | LRIB | LTH1 | LTH2 |
    +--------+------+--------+-------+------+------+------+------+------+
    | THI1   | THN1 |  THI2  |  THN2 | THI3 | THN3 |      |      |      |
    +--------+------+--------+-------+------+------+------+------+------+

    """
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.caero_body_id = np.array([], dtype='int32')

        self.lrsb = np.array([], dtype='int32')
        self.lrib = np.array([], dtype='int32')
        self.lth = np.zeros((0, 2), dtype='int32')
        self.aspect_ratio = np.array([], dtype='float64')
        self.width = np.array([], dtype='float64')
        self.orientation = np.array([], dtype='|U4')
        self.thi = np.zeros((0, 3), dtype='int32')
        self.thn = np.zeros((0, 3), dtype='int32')

    def add(self, pid: int, orient: str, width: float, AR: float,
            thi: list[int], thn: list[int],
            lrsb: Optional[int]=None,
            lrib: Optional[int]=None,
            lth: Optional[int]=None,
            comment: str='') -> int:
        """
        Creates a PAERO2 card, which defines additional cross-sectional
        properties for the CAERO2 geometry.

        Parameters
        ----------
        pid : int
            PAERO1 id
        orient : str
            Orientation flag. Type of motion allowed for bodies. Refers
            to the aerodynamic coordinate system of ACSID. See AERO entry.
            valid_orientations = {Z, Y, ZY}
        width : float
            Reference half-width of body and the width of the constant
            width interference tube
        AR : float
            Aspect ratio of the interference tube (height/width)
        thi / thn : list[int]
            The first (thi) and last (thn) interference element of a body
            to use the theta1/theta2 array
        lrsb : int; default=None
            int : AEFACT id containing a list of slender body half-widths
                  at the end points of the slender body elements
            None : use width
        lrib : int; default=None
            int : AEFACT id containing a list of interference body
                  half-widths at the end points of the interference elements
            None : use width
        lth : list[int, int]; default=None
            AEFACT ids for defining theta arrays for interference calculations
            for theta1/theta2; length=2
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, orient, width, AR, thi, thn,
                           lrsb, lrib, lth, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PAERO2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        orient = string(card, 2, 'orient')
        width = double(card, 3, 'width')
        AR = double(card, 4, 'AR')
        lrsb = integer_or_blank(card, 5, 'lrsb')
        lrib = integer_or_blank(card, 6, 'lrib')
        lth1 = integer_or_blank(card, 7, 'lth1', default=0)
        lth2 = integer_or_blank(card, 8, 'lth2', default=0)
        thi = []
        thn = []
        list_fields = [interpret_value(field, card) for field in card[9:]]
        nfields = len(list_fields)
        lth = [lth1, lth2]
        for i in range(9, 9 + nfields, 2):
            thi.append(integer(card, i, 'lth'))
            thn.append(integer(card, i + 1, 'thn'))
        #return PAERO2(pid, orient, width, AR, thi, thn,
                      #lrsb=lrsb, lrib=lrib, lth=lth,
                      #comment=comment)
        assert len(thi) <= 3, thi
        self.cards.append((pid, orient, width, AR, thi, thn,
                           lrsb, lrib, lth, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        #caero_body_id = np.zeros(ncards, dtype='int32')
        lrsb = np.zeros(ncards, dtype='int32')
        lrib = np.zeros(ncards, dtype='int32')
        lth = np.zeros((ncards, 2), dtype='int32')

        aspect_ratio = np.zeros(ncards, dtype='float64')
        width = np.zeros(ncards, dtype='float64')
        orientation = np.zeros(ncards, dtype='|U4')
        thi = np.zeros((ncards, 3), dtype='int32')
        thn = np.zeros((ncards, 3), dtype='int32')

        for icard, card in enumerate(self.cards):
            (pid, orienti, widthi, AR, thii, thni,
             lrsbi, lribi, lthi, comment) = card
            property_id[icard] = pid
            aspect_ratio[icard] = AR
            width[icard] = widthi
            orientation[icard] = orienti
            thi[icard, :len(thii)] = thii
            thn[icard, :len(thni)] = thni
            lrsbi = lrsbi if lrsbi is not None else 0
            lribi = lribi if lribi is not None else 0
            lthi = lthi if lthi is not None else 0
            lrsb[icard] = lrsbi
            lrib[icard] = lribi
            if isinstance(lthi, integer_types):
                lthi = [lthi]
            lth[icard, :len(lthi)] = lthi
            #x = 1
            #if ncardsi == 0:
                #continue
        self._save(property_id, aspect_ratio, width, orientation,
                   thi, thn, lrsb, lrib, lth)
        self.sort()
        self.cards = []

    def _save(self, property_id, aspect_ratio, width, orientation,
              thi, thn, lrsb, lrib, lth):
        if len(self.property_id):
            property_id = np.hstack([self.property_id, property_id])
            aspect_ratio = np.hstack([self.aspect_ratio, aspect_ratio])
            width = np.hstack([self.width, width])
            orientation = np.hstack([self.orientation, orientation])
            thi = np.vstack([self.thi, thi])
            thn = np.vstack([self.thn, thn])
            lrsb = np.hstack([self.lrsb, lrsb])
            lrib = np.hstack([self.lrib, lrib])
            lth = np.vstack([self.lth, lth])

        self.property_id = property_id
        self.aspect_ratio = aspect_ratio
        self.width = width
        self.orientation = orientation
        self.thi = thi
        self.thn = thn
        self.lrsb = lrsb
        self.lrib = lrib
        self.lth = lth
        self.n = len(property_id)

    def validate(self) -> None:
        if self.n == 0:
            return
        nrows = len(self.property_id)
        assert self.lth.shape == (nrows, 2), self.lth.shape
        assert self.thi.shape == (nrows, 3), self.thi.shape
        assert self.thn.shape == (nrows, 3), self.thn.shape

    def __apply_slice__(self, prop: PAERO2, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.aspect_ratio = self.aspect_ratio[i]
        prop.width = self.width[i]
        prop.orientation = self.orientation[i]
        prop.thi = self.thi[i]
        prop.thn = self.thn[i]
        prop.lrsb = self.lrsb[i]
        prop.lrib = self.lrib[i]
        prop.lth = self.lth[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        #geom_check(
            #missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.element_id, ucaero_ids),
        #)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        paero_ids = array_str(self.property_id, size=size)
        lrsbs = array_str(self.lrsb, size=size)
        lribs = array_str(self.lrib, size=size)

        for pid, orient, width, AR, lrsb, lrib, lth, thi, thn in zip(
            paero_ids, self.orientation, self.width, self.aspect_ratio,
            lrsbs, lribs, self.lth, self.thi, self.thn):
            list_fields = ['PAERO2', pid, orient, width,
                           AR, lrsb, lrib] + list(lth)
            for (thii, thni) in zip(thi, thn):
                list_fields += [thii, thni]
            bdf_file.write(print_card(list_fields))
        return


class PAERO3(PAERO):
    """
    Defines the number of Mach boxes in the flow direction and the
    location of cranks and control surfaces of a Mach box lifting
    surface.

    +--------+------+------+-------+------+-----+------+------+------+
    |    1   |   2  |   3  |   4   |   5  |  6  |   7  |   8  |  9   |
    +========+======+======+=======+======+=====+======+======+======+
    | PAERO3 |  PID | NBOX | NCTRL |      |  X5 |  Y5  |  X6  |  Y6  |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  X7  |  Y7  |   X8  |  Y8  |  X9 |  Y9  |  X10 |  Y10 |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  X11 |  Y11 |  X12  |  Y12 |     |      |      |      |
    +--------+------+------+-------+------+-----+------+------+------+
    | PAERO3 | 2001 |  15  |   1   |      | 0.  |  65. |      |      |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  78. |  65. |  108. |  65. | 82. | 97.5 | 112. | 97.5 |
    +--------+------+------+-------+------+-----+------+------+------+
    |        |  86. | 130. |  116. | 130. |     |      |      |      |
    +--------+------+------+-------+------+-----+------+------+------+

    """
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')
        self.caero_body_id = np.array([], dtype='int32')

    def add(self, pid: int, nbox: int, ncontrol_surfaces: int,
            x: list[float], y: list[float], comment: str='') -> int:
        """
        Creates a PAERO3 card, which defines the number of Mach boxes
        in the flow direction and the location of cranks and control
        surfaces of a Mach box lifting surface.

        Parameters
        ----------
        pid : int
            PAERO1 id
        nbox : int
            Number of Mach boxes in the flow direction; 0 < nbox < 50
        ncontrol_surfaces : int
            Number of control surfaces. (0, 1, or 2)
        x / y : list[float, None]
            float : locations of points 5 through 12, which are in the
            aerodynamic coordinate system, to define the cranks and
            control surface geometry.
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((pid, nbox, ncontrol_surfaces, x, y, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PAERO1 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        """
        Adds a PAERO3 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        nbox = integer(card, 2, 'nbox')
        ncontrol_surfaces = integer(card, 3, 'ncontrol_surfaces')
        x = []
        y = []
        nfields = card.nfields

        j = 5
        for i in range(5, nfields, 2):
            xi = double_or_blank(card, i, 'x%d' % j)
            yi = double_or_blank(card, i + 1, 'y%d' % j)
            x.append(xi)
            y.append(yi)
            j += 1
        #return PAERO3(pid, nbox, ncontrol_surfaces, x, y, comment=comment)
        self.cards.append((pid, nbox, ncontrol_surfaces, x, y, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        nbox = np.zeros(ncards, dtype='int32')
        ncontrol_surface = np.zeros(ncards, dtype='int32')
        x = np.full((ncards, 8), np.nan, dtype='float64')
        y = np.full((ncards, 8), np.nan, dtype='float64')
        for icard, card in enumerate(self.cards):
            (pid, nboxi, ncontrol_surfacei, xi, yi, comment) = card
            property_id[icard] = pid
            nbox[icard] = nboxi
            ncontrol_surface[icard] = ncontrol_surfacei
            x[icard, :len(xi)] = xi
            y[icard, :len(yi)] = yi
        self._save(property_id, nbox, ncontrol_surface, x, y)
        self.sort()
        self.cards = []

    def _save(self, property_id, nbox, ncontrol_surface, x, y):
        assert len(self.property_id) == 0, self.property_id
        self.property_id = property_id
        self.nbox = nbox
        self.ncontrol_surface = ncontrol_surface

        #[x5, x6, x7, x8], [x9, x10, x11, x12]
        self.x = x
        self.y = y
        self.n = len(property_id)

    def __apply_slice__(self, prop: PAERO3, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        prop.nbox = self.nbox[i]
        prop.ncontrol_surface = self.ncontrol_surface[i]
        prop.x = self.x[i]
        prop.y = self.y[i]

    #def __apply_slice__(self, prop: PAERO3, i: np.ndarray) -> None:
        #prop.n = len(i)
        #prop.property_id = self.property_id[i]
        #prop.caero_body_id = self.caero_body_id[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        #geom_check(
            #missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.element_id, ucaero_ids),
        #)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        paero_ids = array_str(self.property_id, size=size)
        nboxes = array_str(self.nbox, size=size)
        ncontrol_surfaces = array_str(self.ncontrol_surface, size=size)

        x = array_float_nan(self.x, size=size, is_double=False)
        y = array_float_nan(self.y, size=size, is_double=False)
        for pid, nboxi, ncontrol_surfacei, xi, yi in zip(paero_ids, nboxes, ncontrol_surfaces, x, y):
            list_fields = ['PAERO3', pid, nboxi, ncontrol_surfacei, None]
            for (xii, yii) in zip(xi, yi):
                list_fields += [xii, yii]
            bdf_file.write(print_card(list_fields))
        return


class PAERO4(PAERO):
    """
    Defines properties of each strip element for Strip theory.

    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |    1   |   2  |   3   |   4    |   5   |   6   |    7   |   8    |    9   |
    +========+======+=======+========+=======+=======+========+========+========+
    | PAERO4 | PID  | CLA   |  LCLA  |  CIRC | LCIRC |  DOC1  |  CAOC1 | GAPOC1 |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |        | DOC2 | CAOC2 | GAPOC2 |  DOC3 | CAOC3 | GAPOC3 |  etc.  |        |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    | PAERO4 | 6001 |   1   |   501  |   0   |   0   |   0.0  |   0.0  |   0.0  |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    |        | 0.50 |  0.25 |  0.02  |  0.53 |  0.24 |   0.0  |        |        |
    +--------+------+-------+--------+-------+-------+--------+--------+--------+
    ## TODO: what happens for DOC4?

    """
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def add(self, pid: int,
            docs: list[float], caocs: list[float], gapocs: list[float],
            cla: int=0, lcla: int=0,
            circ: int=0, lcirc: int=0, comment: str='') -> int:
        """
        Parameters
        ----------
        PID : int
            Property identification number. (Integer > 0)
        CLA : int; default=0
            Select Prandtl-Glauert correction. (Integer = -1, 0, 1)
            -1 Compressibility correction made to lift curve slope data for a reference Mach number.
            0  No correction and no list needed. (Default)
            +1 No correction and lift curve slope provided by a list as a
               function of strip location and Mach number.
        LCLA : int
            ID number of the AEFACT entry that lists the lift curve slope
            on all strips for each Mach number on the MKAEROi entry. See
            Remark 2 below. (Integer = 0 if CLA = 0, > 0 if CLA ≠ 0)
        CIRC : int; default=0
            Select Theodorsen’s function C(k) or the number of exponential
            coefficients used to approximate C(k).
            (Integer = 0, 1, 2, 3; Must be zero if CLA ≠ 0.)
            0 Theodorsen function.
            1, 2, 3 Approximate function with b0, b1, β1, ..., bn, βn n = 1, 2, 3.
        LCIRC : int
            Identification number of the AEFACT entry that lists the b, β values
            for each Mach number. See Remark 3, 4, and 5 below; variable b’s
            and β’s for each mi on the MKAEROi entry.
            (Integer = 0 if CIRC = 0, > 0 if CIRC ≠ 0)
        DOCi : list[float]
            d/c = distance of the control surface hinge aft of the quarter-chord
            divided by the strip chord (Real ≥ 0.0)
        CAOCi : list[float]
            ca/c = control surface chord divided by the strip chord. (Real ≥ 0.0)
        GAPOCi : list[float]
            g/c = control surface gap divided by the strip chord. (Real ≥ 0.0)

        """
        card = (pid, docs, caocs, gapocs,
                cla, lcla, circ, lcirc, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PAERO4 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'pid')
        cla = integer_or_blank(card, 2, 'cla', default=0)
        lcla = integer_or_blank(card, 3, 'lcla', default=0) # ???

        circ = integer_or_blank(card, 4, 'circ', default=0)
        lcirc = integer_or_blank(card, 5, 'lcirc', default=0) # ???
        nfields = card.nfields

        j = 0
        docs = []
        caocs = []
        gapocs = []
        for i in range(6, nfields, 3):
            doc = double(card, i, 'doc_%d' % j)
            caoc = double(card, i + 1, 'caoc_%d' % j)
            gapoc = double(card, i + 2, 'gapoc_%d' % j)
            docs.append(doc)
            caocs.append(caoc)
            gapocs.append(gapoc)
            j += 1
        #return PAERO4(pid, docs, caocs, gapocs,
                      #cla=cla, lcla=lcla, circ=circ, lcirc=lcirc, comment=comment)
        card = (pid, docs, caocs, gapocs,
                cla, lcla, circ, lcirc, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        #caero_body_id = np.zeros(ncards, dtype='int32')
        #docs = np.zeros(ncards, dtype='float64')
        #caocs = np.zeros(ncards, dtype='float64')
        #gapocs = np.zeros(ncards, dtype='float64')
        docs = []
        caocs = []
        gapocs = []
        cla = np.zeros(ncards, dtype='int32')
        lcla = np.zeros(ncards, dtype='int32')
        circ = np.zeros(ncards, dtype='int32')
        lcirc = np.zeros(ncards, dtype='int32')
        ndoc = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (pid, docsi, caocsi, gapocsi,
             clai, lclai, circi, lcirci, comment) = card
            property_id[icard] = pid

            ndoci = len(docs)
            docs.append(docsi)
            caocs.append(caocsi)
            gapocs.append(gapocsi)
            ndoc[icard] = ndoci

            cla[icard] = clai
            lcla[icard] = lclai
            circ[icard] = circi
            lcirc[icard] = lcirci
        self._save(property_id, ndoc, docs, caocs, gapocs,
                   cla, lcla, circ, lcirc)
        self.sort()
        self.cards = []

    def _save(self, property_id, ndoc, docs, caocs, gapocs,
              cla, lcla, circ, lcirc):
        assert len(self.property_id) == 0, self.property_id
        self.property_id = property_id
        self.ndoc = ndoc
        self.docs = docs
        self.caocs = caocs
        self.gapocs = gapocs
        self.cla = cla
        self.lcla = lcla
        self.circ = circ
        self.lcirc = lcirc
        self.n = len(property_id)

    @property
    def idoc(self) -> np.ndarray:
        return make_idim(self.n, self.ndoc)

    def __apply_slice__(self, prop: PAERO4, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]

        idoc = self.idoc
        prop.docs = hslice_by_idim(i, idoc, self.docs)
        prop.caocs = hslice_by_idim(i, idoc, self.caocs)
        prop.gapocs = hslice_by_idim(i, idoc, self.gapocs)
        prop.ndoc = self.ndoc[i]

        prop.cla = self.cla[i]
        prop.lcla = self.lcla[i]
        prop.circ = self.circ[i]
        prop.lcirc = self.lcirc[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        #geom_check(
            #missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.element_id, ucaero_ids),
        #)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        paero_ids = array_str(self.property_id, size=size)
        cla_ = array_str(self.cla, size=size)
        lcla_ = array_str(self.lcla, size=size)
        circ_ = array_str(self.circ, size=size)
        lcirc_ = array_str(self.lcirc, size=size)

        for pid, cla, lcla, circ, lcirc, (idoc1, idoc2) in zip(
                paero_ids, cla_, lcla_, circ_, lcirc_, self.idoc):
            list_fields = ['PAERO4', pid, cla, lcla, circ, lcirc]
            docs = self.docs[idoc1:idoc2]
            caocs = self.caocs[idoc1:idoc2]
            gapocs = self.gapocs[idoc1:idoc2]
            for doc, caoc, gapoc in zip(docs, caocs, gapocs):
                list_fields += [doc, caoc, gapoc]
            bdf_file.write(print_card(list_fields))
        return


class PAERO5(PAERO):
    """
    +--------+-------+--------+--------+---------+-------+-------+-------+
    |   1    |   2   |    3   |   4    |    5    |   6   |   7   |   8   |
    +========+=======+========+========+=========+=======+=======+=======+
    | PAERO5 |  PID  | NALPHA | LALPHA |  NXIS   | LXIS  | NTAUS | LTAUS |
    +--------+-------+--------+--------+---------+-------+-------+-------+
    |        | CAOC1 | CAOC2  | CAOC3  |  CAOC4  | CAOC5 |       |       |
    +--------+-------+--------+--------+---------+-------+-------+-------+
    | PAERO5 | 7001  |   1    |  702   |    1    | 701   |   1   |  700  |
    +--------+-------+--------+--------+---------+-------+-------+-------+
    |        |  0.0  |  0.0   |  5.25  | 3.99375 |  0.0  |       |       |
    +--------+-------+--------+--------+---------+-------+-------+-------+
    """
    _id_name = 'property_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.property_id = np.array([], dtype='int32')

    def add(self, pid: int, caoci: list[float],
            nalpha: int=0, lalpha: int=0,
            nxis: int=0, lxis: int=0,
            ntaus: int=0, ltaus: int=0,
            comment='') -> int:
        """Creates a PAERO5 card"""
        card = (pid, caoci,
                nalpha, lalpha, nxis, lxis,
                ntaus, ltaus,
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a PAERO5 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        pid = integer(card, 1, 'property_id')
        nalpha = integer_or_blank(card, 2, 'nalpha', default=0)
        lalpha = integer_or_blank(card, 3, 'lalpha', default=0)

        nxis = integer_or_blank(card, 4, 'nxis', default=0)
        lxis = integer_or_blank(card, 5, 'lxis', default=0)

        ntaus = integer_or_blank(card, 6, 'ntaus', default=0)
        ltaus = integer_or_blank(card, 7, 'ltaus', default=0)

        caoci = []
        for n, i in enumerate(range(9, len(card))):
            ca = double(card, i, 'ca/ci_%i' % (n+1))
            caoci.append(ca)
        #return PAERO5(pid, caoci,
                      #nalpha=nalpha, lalpha=lalpha, nxis=nxis, lxis=lxis,
                      #ntaus=ntaus, ltaus=ltaus,
                      #comment=comment)
        card = (pid, caoci,
                nalpha, lalpha, nxis, lxis,
                ntaus, ltaus,
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        property_id = np.zeros(ncards, dtype='int32')
        ncaoci = np.zeros(ncards, dtype='int32')
        caoci = []
        nalpha = np.zeros(ncards, dtype='int32')
        lalpha = np.zeros(ncards, dtype='int32')
        nxis = np.zeros(ncards, dtype='int32')
        lxis = np.zeros(ncards, dtype='int32')
        ntaus = np.zeros(ncards, dtype='int32')
        ltaus = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (pid, caocii,
             nalphai, lalphai, nxisi, lxisi,
             ntausi, ltausi,
             comment) = card
            property_id[icard] = pid

            #ndoci = len(docs)
            caoci.append(caocii)
            ncaoci[icard] = len(caocii)
            #ndoc[icard] = ndoci

            nalpha[icard] = nalphai
            lalpha[icard] = lalphai
            nxis[icard] = nxisi
            lxis[icard] = lxisi
            ntaus[icard] = ntausi
            ltaus[icard] = ltausi

        caoci = np.array(caoci, dtype='float64')
        self._save(property_id, ncaoci, caoci, nalpha, lalpha, nxis, lxis, ntaus, ltaus)
        self.sort()
        self.cards = []

    def _save(self, property_id, ncaoci, caoci, nalpha, lalpha,
              nxis, lxis, ntaus, ltaus):
        assert len(self.property_id) == 0, self.property_id
        self.property_id = property_id
        #self.ndoc = ndoc
        self.ncaoci = ncaoci
        self.caoci = caoci

        self.nalpha = nalpha
        self.lalpha = lalpha
        self.nxis = nxis
        self.lxis = lxis
        self.ntaus = ntaus
        self.ltaus = ltaus
        self.n = len(property_id)

    def __apply_slice__(self, prop: PAERO5, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.property_id = self.property_id[i]
        #prop.docs = self.docs[i]
        #prop.caocs = self.caocs[i]
        #prop.gapocs = self.gapocs[i]

        icaoci = self.icaoci
        prop.docs = hslice_by_idim(i, icaoci, self.docs)
        prop.caocs = hslice_by_idim(i, icaoci, self.caocs)
        prop.gapocs = hslice_by_idim(i, icaoci, self.gapocs)
        prop.ndoc = self.ndoc[i]

        prop.cla = self.cla[i]
        prop.lcla = self.lcla[i]
        prop.circ = self.circ[i]
        prop.lcirc = self.lcirc[i]

    @property
    def icaoci(self) -> np.ndarray:
        return make_idim(self.n, self.ncaoci)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        #ucaero_ids = np.unique(self.caero_id)

        #set1_ids = np.unique(set1_ids)
        #geom_check(
            #missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.element_id, ucaero_ids),
        #)

    @parse_property_check
    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        print_card = get_print_card_8_16(size)

        paero_ids = array_str(self.property_id, size=size)
        lalpha_ = array_str(self.lalpha, size=size)
        nalpha_ = array_str(self.nalpha, size=size)

        nxis_ = array_str(self.nxis, size=size)
        lxis_ = array_str(self.lxis, size=size)
        ntaus_ = array_str(self.ntaus, size=size)
        ltaus_ = array_str(self.ltaus, size=size)

        for pid, nalpha, lalpha, nxis, lxis, ntau, ltau, (caoci0, caoci1) in zip(
                paero_ids, nalpha_, lalpha_,
                nxis_, lxis_,
                ntaus_, ltaus_, self.icaoci):
            caoci = self.caoci[caoci0:caoci1].ravel()
            list_fields = ['PAERO5', pid, nalpha, lalpha, nxis,
                           lxis, ntau, ltau] + list(caoci)
            bdf_file.write(print_card(list_fields))
        return


class AELIST(VectorizedBaseCard):
    """
    Defines a list of aerodynamic elements to undergo the motion prescribed
    with the AESURF Bulk Data entry for static aeroelasticity.

    +---------+------+------+------+------+------+------+------+------+
    |    1    |   2  |   3  |  4   |   5  |   6  |  7   |  8   |  9   |
    +=========+======+======+======+======+======+======+======+======+
    |  AELIST |  SID |  E1  |  E2  |  E3  |  E4  |  E5  |  E6  |  E7  |
    +---------+------+------+------+------+------+------+------+------+
    |         |  E8  | etc. |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+
    |  AELIST |  75  | 1001 | THRU | 1075 | 1101 | THRU | 1109 | 1201 |
    +---------+------+------+------+------+------+------+------+------+
    |         | 1202 |      |      |      |      |      |      |      |
    +---------+------+------+------+------+------+------+------+------+

    Notes
    -----
    1. These entries are referenced by the AESURF entry.
    2. When the THRU option is used, all intermediate grid points must exist.
       The word THRU may not appear in field 3 or 9 (2 or 9 for continuations).
    3. Intervening blank fields are not allowed.

    """
    _id_name = 'aelist_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aelist_id = np.array([], dtype='int32')
        self.nelements = np.array([], dtype='int32')
        self.elements = np.array([], dtype='int32')

    def add(self, sid: int, elements: list[int], comment: str='') -> int:
        """
        Creates an AELIST card, which defines the aero boxes for
        an AESURF/SPLINEx.

        Parameters
        ----------
        sid : int
            unique id
        elements : list[int, ..., int]
            list of box ids
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, elements, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AELIST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        elements = fields(integer_or_string, card, 'eid', i=2, j=len(card))
        #print(f'sid={sid} elements={elements}')
        self.cards.append((sid, elements, comment))
        self.n += 1
        #return AELIST(sid, elements, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aelist_id = np.zeros(ncards, dtype='int32')
        nelements = np.zeros(ncards, dtype='int32')
        #elements = np.array([], dtype='int32')

        all_elements = []
        for icard, card in enumerate(self.cards):
            sid, elementsi, comment = card
            aelist_id[icard] = sid
            elements2 = expand_thru(elementsi)
            #elements2.sort()
            nelements[icard] = len(elements2)
            all_elements.extend(elements2)
        elements = np.array(all_elements, dtype='int32')
        self._save(aelist_id, nelements, elements)
        self.clean_ids()
        self.sort()
        self.cards = []

    def clean_ids(self):
        elements = []
        nelements = np.zeros(self.nelements.shape, dtype=self.nelements.dtype)
        for ie, (ielement0, ielement1) in enumerate(self.ielement):
            element = self.elements[ielement0:ielement1]
            uelement = np.unique(element)
            nelement = len(uelement)
            elements.append(uelement)
            nelements[ie] = nelement
        self.elements = np.hstack(elements)
        self.nelements = nelements

    def _save(self, aelist_id, nelements, elements):
        if len(self.aelist_id):
            aelist_id = np.hstack([self.aelist_id, aelist_id])
            elements = np.hstack([self.elements, elements])
            nelements = np.hstack([self.nelements, nelements])
        self.aelist_id = aelist_id
        self.elements = elements
        self.nelements = nelements
        self.n = len(aelist_id)

    def __apply_slice__(self, card: AELIST, i: np.ndarray) -> None:
        card.n = len(i)
        card.aelist_id = self.aelist_id[i]

        ielement = self.ielement
        card.elements = hslice_by_idim(i, ielement, self.elements)
        card.nelements = self.nelements[i]

    def slice_card_by_aelist_id(self, aelist_id: np.ndarray) -> AELIST:
        assert self.n > 0, self.n
        assert len(self.aelist_id) > 0, self.aelist_id
        i = self.index(aelist_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    @property
    def ielement(self) -> np.ndarray:
        return make_idim(self.n, self.nelements)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aelist_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        elements_ = array_str(self.elements, size=size).tolist()
        for sid, (ielement0, ielement1) in zip(self.aelist_id, self.ielement):
            elements = elements_[ielement0:ielement1]
            elementsi = self.elements[ielement0:ielement1]
            elementsi.sort()
            nelement = ielement1 - ielement0
            if elementsi[0] + nelement == elementsi[-1] - 1:
                list_fields = ['AELIST', sid, elementsi[0], 'THRU', elementsi[-1] - 1]
            else:
                list_fields = ['AELIST', sid] + elements
            bdf_file.write(print_card(list_fields))
        return


class AELINK(VectorizedBaseCard):
    r"""
    Defines relationships between or among AESTAT and AESURF entries, such
    that:

    .. math:: u^D + \Sigma_{i=1}^n C_i u_i^I = 0.0

    +--------+-------+-------+--------+------+-------+----+-------+----+
    |   1    |   2   |   3   |   4    |   5  |   6   |  7 |   8   |  9 |
    +========+=======+=======+========+======+=======+====+=======+====+
    | AELINK |  ID   | LABLD | LABL1  |  C1  | LABL2 | C2 | LABL3 | C3 |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    |        | LABL4 |  C4   |  etc.  |      |       |    |       |    |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    | AELINK |  10   | INBDA |  OTBDA | -2.0 |       |    |       |    |
    +--------+-------+-------+--------+------+-------+----+-------+----+
    """
    _id_name = 'aelink_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aelink_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.independent_labels = np.array([], dtype='|U8')
        self.linking_coefficients = np.array([], dtype='float64')
        self.nindependent_labels = np.array([], dtype='int32')

    def add(self, aelink_id: int, label: str,
            independent_labels: list[str],
            linking_coefficients: list[float],
            comment: str='') -> int:
        """
        Creates an AELINK card, which defines an equation linking
        AESTAT and AESURF cards

        Parameters
        ----------
        aelink_id : int
            unique id
        label : str
            name of the dependent AESURF card
        independent_labels : list[str, ..., str]
            name for the independent variables (AESTATs)
        linking_coefficients : list[float]
            linking coefficients
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((aelink_id, label, independent_labels,
                           linking_coefficients, comment))
        assert len(independent_labels) == len(linking_coefficients)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AELINK card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aelink_id = integer_or_string(card, 1, 'ID')
        label = string(card, 2, 'label')
        independent_labels = []
        linking_coefficients = []

        list_fields = [interpret_value(field, card) for field in card[3:]]
        assert len(list_fields) % 2 == 0, 'list_fields=%s' % list_fields
        for i in range(0, len(list_fields), 2):
            independent_label = list_fields[i]
            linking_coefficent = list_fields[i + 1]
            independent_labels.append(independent_label)
            linking_coefficients.append(linking_coefficent)
        #return AELINK(aelink_id, label, independent_labels, linking_coefficients,
                      #comment=comment)
        assert len(independent_labels) == len(linking_coefficients)
        self.cards.append((aelink_id, label, independent_labels,
                           linking_coefficients, comment))
        self.n += 1
        #return AELIST(sid, elements, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aelink_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        independent_labels = []
        linking_coefficients = []
        nindependent_labels = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (aelink_idi, labeli, independent_labelsi,
             linking_coefficientsi, comment) = card
            if aelink_idi == 'ALWAYS':
                aelink_idi = 0
            aelink_id[icard] = aelink_idi
            label[icard] = labeli
            independent_labels.extend(independent_labelsi)
            linking_coefficients.extend(linking_coefficientsi)
            assert len(independent_labelsi) == len(linking_coefficientsi)
            nindependent_labels[icard] = len(independent_labelsi)
            #elements2.sort()
        independent_labels = np.array(independent_labels, dtype='|U8')
        linking_coefficients = np.array(linking_coefficients, dtype='float64')
        self._save(aelink_id, label,
                   nindependent_labels, independent_labels, linking_coefficients)
        self.sort()
        self.cards = []

    def _save(self, aelink_id, label,
              nindependent_labels, independent_labels, linking_coefficients):
        if len(self.aelink_id):
            aelink_id = np.hstack([self.aelink_id, aelink_id])
            label = np.hstack([self.label, label])
            independent_labels = np.hstack([self.independent_labels, independent_labels])
            nindependent_labels = np.hstack([self.nindependent_labels, nindependent_labels])
            linking_coefficients = np.hstack([self.linking_coefficients, linking_coefficients])

        self.aelink_id = aelink_id
        self.label = label
        self.independent_labels = independent_labels
        self.nindependent_labels = nindependent_labels
        self.linking_coefficients = linking_coefficients
        self.n = len(aelink_id)

    @property
    def icoeff(self) -> np.ndarray:
        return make_idim(self.n, self.nindependent_labels)

    def __apply_slice__(self, card: AELINK, i: np.ndarray) -> None:
        card.n = len(i)
        card.aelink_id = self.aelink_id[i]

        icoeff = self.icoeff
        card.independent_labels = hslice_by_idim(i, icoeff, self.independent_labels)
        card.linking_coefficients = hslice_by_idim(i, icoeff, self.linking_coefficients)
        card.nindependent_labels = self.nindependent_labels[i]

    def slice_card_by_aelink_id(self, aelink_id: np.ndarray) -> AELINK:
        assert self.n > 0, self.n
        assert len(self.aelink_id) > 0, self.aelink_id
        i = self.index(aelink_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aelink_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        aelink_ids = array_str(self.aelink_id, size=size)
        for aelink_id, label, (icoeff0, icoeff1) in zip(aelink_ids, self.label, self.icoeff):

            independent_labels = self.independent_labels[icoeff0:icoeff1]
            linking_coefficients = self.linking_coefficients[icoeff0:icoeff1]

            list_fields = ['AELINK', aelink_id, label]
            for (ivar, ival) in zip(independent_labels, linking_coefficients):
                list_fields += [ivar, ival]
            bdf_file.write(print_card(list_fields))
        return


class AEFACT(VectorizedBaseCard):
    """
    Defines real numbers for aeroelastic analysis.

    +--------+-----+----+--------+-----+----+----+----+----+
    |   1    |   2 |  3 |    4   |  5  |  6 |  7 | 8  |  9 |
    +========+=====+====+========+=====+====+====+====+====+
    | AEFACT | SID | D1 |   D2   | D3  | D4 | D5 | D6 | D7 |
    +--------+-----+----+--------+-----+----+----+----+----+
    |        | D8  | D9 |  etc.  |     |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+
    | AEFACT | 97  |.3  |  0.7   | 1.0 |    |    |    |    |
    +--------+-----+----+--------+-----+----+----+----+----+

    TODO: Are these defined in percentages and thus,
          should they be normalized if they are not?

    """
    _id_name = 'aefact_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aefact_id = np.array([], dtype='int32')
        self.nfractions = np.array([], dtype='int32')
        self.fractions = np.array([], dtype='float64')

    #def __len__(self):
        #return len(self.aefact_id)

    def slice_card_by_aefact_id(self, aefact_id: np.ndarray) -> AEFACT:
        assert self.n > 0, self.n
        assert len(self.aefact_id) > 0, self.aefact_id
        i = self.index(aefact_id)
        #cls_obj = cls(self.model)
        #cls_obj.__apply_slice__(self, i)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def add(self, sid: int, fractions: list[float], comment: str='') -> int:
        """
        Creates an AEFACT card, which is used by the CAEROx / PAEROx card
        to adjust the spacing of the sub-paneleing (and grid point
        paneling in the case of the CAERO3).

        Parameters
        ----------
        sid : int
            unique id
        fractions : list[float, ..., float]
            list of percentages
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, fractions, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> None:
        sid = integer(card, 1, 'sid')
        fractions = []
        for i in range(2, len(card)):
            fraction = double(card, i, 'factor_%d' % (i - 1))
            fractions.append(fraction)
        assert len(card) > 2, 'len(AEFACT card) = %i\n%s' % (len(card), card)
        self.cards.append((sid, fractions, comment))
        self.n += 1
        #return AEFACT(sid, fractions, comment=comment)
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aefact_id = np.zeros(ncards, dtype='int32')
        nfractions = np.zeros(ncards, dtype='int32')
        #self.fractions = np.array([], dtype='int32')

        all_fractions = []
        for icard, card in enumerate(self.cards):
            sid, fractions, comment = card
            aefact_id[icard] = sid
            nfractions[icard] = len(fractions)
            all_fractions.extend(fractions)
        fractions = np.array(all_fractions, dtype='float64')
        self._save(aefact_id, nfractions, fractions)
        self.cards = []

    def _save(self, aefact_id, nfractions, fractions) -> None:
        if len(self.aefact_id):
            aefact_id = np.hstack([self.aefact_id, aefact_id])
            nfractions = np.hstack([self.nfractions, nfractions])
            fractions = np.hstack([self.fractions, fractions])
        self.aefact_id = aefact_id
        self.nfractions = nfractions
        self.fractions = fractions

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    @property
    def ifraction(self) -> np.ndarray:
        return make_idim(self.n, self.nfractions)

    def __apply_slice__(self, prop: AEFACT, i: np.ndarray) -> None:
        prop.n = len(i)
        prop.aefact_id = self.aefact_id[i]

        ifraction = self.ifraction
        prop.fractions = hslice_by_idim(i, ifraction, self.fractions)
        prop.nfractions = self.nfractions[i]

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aefact_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        for sid, (ifraction0, ifraction1) in zip(self.aefact_id, self.ifraction):
            factors = self.fractions[ifraction0:ifraction1].tolist()
            list_fields = ['AEFACT', sid] + factors
            bdf_file.write(print_card(list_fields))
        return


class FLFACT(VectorizedBaseCard):
    """
    +--------+-----+----+------+-----+----+----+----+----+
    |   1    |  2  |  3 |   4  |  5  | 6  | 7  | 8  | 9  |
    +========+=====+====+======+=====+====+====+====+====+
    | FLFACT | SID | F1 | F2   | F3  | F4 | F5 | F6 | F7 |
    +--------+-----+----+------+-----+----+----+----+----+
    |        | F8  | F9 | etc. |     |    |    |    |    |
    +--------+-----+----+------+-----+----+----+----+----+
    | FLFACT | 97  | .3 |  .7  | 3.5 |    |    |    |    |
    +--------+-----+----+------+-----+----+----+----+----+

    # delta quantity approach

    +--------+-----+-------+------+-------+----+--------+
    |   1    |  2  |  3    |   4  |   5   | 6  |   7    |
    +========+=====+=======+======+=======+====+========+
    | FLFACT | SID | F1    | THRU | FNF   | NF |  FMID  |
    +--------+-----+-------+------+-------+----+--------+
    | FLFACT | 201 | 0.200 | THRU | 0.100 | 11 | 0.1333 |
    +--------+-----+-------+------+-------+----+--------+
    """
    _id_name = 'flfact_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.flfact_id = np.array([], dtype='int32')
        self.nfactors = np.array([], dtype='int32')
        self.factors = np.array([], dtype='float64')

    def slice_card_by_flfact_id(self, flfact_id: np.ndarray) -> FLFACT:
        assert self.n > 0, self.n
        assert len(self.flfact_id) > 0, self.flfact_id
        i = self.index(flfact_id)
        cls_obj = self.slice_card_by_index(i)
        assert cls_obj.n > 0, cls_obj
        return cls_obj

    def add(self, sid: int, factors: list[float], comment: str='') -> int:
        """
        Creates an FLFACT card, which defines factors used for flutter
        analysis.  These factors define either:
         - density
         - mach
         - velocity
         - reduced frequency
        depending on the FLUTTER method chosen (e.g., PK, PKNL, PKNLS)

        Parameters
        ----------
        sid : int
            the id of a density, reduced_frequency, mach, or velocity table
            the FLUTTER card defines the meaning
        factors : varies
            values : list[float, ..., float]
                list of factors
            list[f1, THRU, fnf, nf, fmid]
                f1 : float
                    first value
                THRU : str
                    the word THRU
                fnf : float
                    second value
                nf : int
                    number of values
                fmid : float; default=(f1 + fnf) / 2.
                    the mid point to bias the array
                TODO: does f1 need be be greater than f2/fnf???
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, factors, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        sid = integer(card, 1, 'sid')
        assert len(card) > 2, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        field3 = double_string_or_blank(card, 3, 'THRU')
        if field3 is None:
            f1 = double(card, 2, 'f1')
            factors = [f1]
            assert len(card) == 3, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
        elif isinstance(field3, float):
            factors = fields(double, card, 'factors', i=2, j=len(card))
        elif isinstance(field3, str) and field3 == 'THRU':
            f1 = double(card, 2, 'f1')
            fnf = double(card, 4, 'fnf')
            nf = integer(card, 5, 'nf')
            fmid_default = (f1 + fnf) / 2.
            fmid = double_or_blank(card, 6, 'fmid', fmid_default)
            assert len(card) <= 7, 'len(FLFACT card)=%s; card=%s' % (len(card), card)
            factors = [f1, 'THRU', fnf, nf, fmid]
        else:
            raise SyntaxError('expected a float or string for FLFACT field 3; value=%r' % field3)
        #return FLFACT(sid, factors, comment=comment)
        self.cards.append((sid, factors, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        flfact_id = np.zeros(ncards, dtype='int32')
        nfactors = np.zeros(ncards, dtype='int32')
        #self.elements = np.array([], dtype='int32')

        all_factors = []
        for icard, card in enumerate(self.cards):
            sid, factors, comment = card
            factors = expand_thru(factors, set_fields=False, sort_fields=False)
            flfact_id[icard] = sid
            nfactors[icard] = len(factors)
            all_factors.extend(factors)
        factors = np.array(all_factors, dtype='float64')
        self._save(flfact_id, nfactors, factors)
        self.cards = []

    def _save(self, flfact_id, nfactors, factors) -> None:
        if len(self.flfact_id):
            flfact_id = np.hstack([self.flfact_id, flfact_id])
            factors = np.hstack([self.factors, factors])
            nfactors = np.hstack([self.nfactors, nfactors])
        self.flfact_id = flfact_id
        self.nfactors = nfactors
        self.factors = factors
        self.n = len(flfact_id)

    def sort(self) -> None:
        uflfact = np.unique(self.flfact_id)
        if np.array_equal(uflfact, self.flfact_id):
            return
        i = np.argsort(self.flfact_id)
        self.__apply_slice__(self, i)

    def __apply_slice__(self, flfact: FLFACT, i: np.ndarray) -> None:
        flfact.flfact_id = self.flfact_id[i]
        ifactor = self.ifactor
        flfact.factors = hslice_by_idim(i, ifactor, self.factors)
        flfact.nfactors = self.nfactors[i]
        flfact.n = len(i)

    @property
    def ifactor(self) -> np.ndarray:
        return make_idim(self.n, self.nfactors)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.flfact_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        for sid, (ifactor0, ifactor1) in zip(self.flfact_id, self.ifactor):
            factors = self.factors[ifactor0:ifactor1].tolist()
            list_fields = ['FLFACT', sid] + factors
            bdf_file.write(print_card(list_fields))
        return


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
            usage: str='BOTH', nelements: int=10,
            melements: int=10, comment: str='') -> int:
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
                nelements, melements, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
                nelements, melements, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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
        for icard, card in enumerate(self.cards):
            (eid, caero, box1, box2, setg, dzi, methodi, usagei,
                nelementi, melementi, comment) = card
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
                   nelement, melement)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, set_id, dz,
              method, usage,
              nelement, melement):
        assert len(self.spline_id) == 0, self.spline_id
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

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spline_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxs = array_str(self.box_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        nelements = array_default_int(self.nelement, default=0, size=size)
        melements = array_default_int(self.melement, default=0, size=size)
        dzs = array_default_float(self.dz, default=0., size=size, is_double=False)
        for eid, caero, (box1, box2), setg, dz, \
            method, usage, nelement, melement in zip(\
                spline_ids, caero_ids, boxs, set_ids, dzs,
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
                #nelements, melements, comment)
        #self.cards.append(card)
        #self.n += 1

    def add(self, eid: int, caero: int,
            box1: int, box2: int, setg: int,
            dz: float=0.0, dtor: float=1.0,
            cid: int=0,
            dthx: float=0.0, dthy: float=0.0,
            usage: str='BOTH',
            comment: str='') -> int:
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
                dthx, dthy, usage, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
                dthx, dthy, usage, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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
        for icard, card in enumerate(self.cards):
            (eid, caero, box1, box2, setg, dzi, dtori, cid,
                dthxi, dthyi, usagei, comment) = card
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
                   coord_id, dtor, dthx, dthy, usage)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, set_id, dz,
                   coord_id, dtor, dthx, dthy, usage) -> None:
        if len(self.spline_id):
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
        all_caero_ids = hstack_lists([model.caero1.element_id, model.caero2.element_id])
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            caero=(all_caero_ids, ucaero_ids),
        )

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spline_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxs = array_str(self.box_id, size=size)
        set_ids = array_str(self.set_id, size=size)

        coord_ids = array_default_int(self.coord_id, default=0, size=size)
        dzs = array_default_float(self.dz, default=0.0, size=size, is_double=False)
        dtors = array_default_float(self.dtor, default=1.0, size=size, is_double=False)
        dthxs = array_default_float(self.dthx, default=0.0, size=size, is_double=False)
        dthys = array_default_float(self.dthy, default=0.0, size=size, is_double=False)
        for eid, caero, (box1, box2), setg, dz, dtor, \
            cid, usage, dthx, dthy in zip(\
                spline_ids, caero_ids, boxs, set_ids, dzs, dtors,
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
            usage: str='BOTH', comment: str='') -> int:
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
                coeffs, usage, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
                coeffs, usage, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        spline_id = np.zeros(ncards, dtype='int32')
        caero_id = np.zeros(ncards, dtype='int32')
        components = np.zeros(ncards, dtype='int32')
        box_id = np.zeros(ncards, dtype='int32')

        nnode = np.zeros(ncards, dtype='int32')
        #ndisp = np.zeros(ncards, dtype='int32')
        #dz = np.zeros(ncards, dtype='float64')

        usage = np.zeros(ncards, dtype='|U4')

        nodes = []
        displacement_components = []
        coeffs = []
        for icard, card in enumerate(self.cards):
            (eid, caero, box_idi, componentsi,
                nodesi, displacement_componentsi,
                coeffsi, usagei, comment) = card
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
                coeffs, usage)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, box_id, components,
                nnode, nodes, displacement_components, coeffs, usage):
        if len(self.spline_id):
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

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spline_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        caero_ids = array_str(self.caero_id, size=size)
        spline_ids = array_str(self.spline_id, size=size)
        boxs = array_str(self.box_id, size=size)
        components = array_str(self.components, size=size)
        nodes = array_str(self.nodes, size=size)
        disp = array_str(self.displacement_components, size=size)

        #spline_id, caero_id, box_id, components,
        #nnodes, nodes, displacement_components, coeffs
        for eid, caero, box_id, component, (inode0, inode1), usage in zip(
                spline_ids, caero_ids, boxs, components,
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
            comment: str='') -> int:
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
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
                comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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

        for icard, card in enumerate(self.cards):
            (eid, caero, aelist, setg, dzi, methodi, usagei,
                nelementi, melementi, ftypei, rcorei,
                comment) = card
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
                   nelement, melement, ftype, rcore)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, aelist_id, set_id, dz,
              method, usage,
              nelement, melement, ftype, rcore):
        assert len(self.spline_id) == 0, self.spline_id
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
        elem.spline_id = self.spline_id[i]
        elem.caero_id = self.caero_id[i]
        elem.aelist_id = self.aelist_id[i]
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

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spline_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

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
            ftype: str='WF2', rcore=None, comment: str='') -> int:
        """Creates a SPLINE5 card"""
        card = (eid, caero, aelist, setg, thx, thy, dz, dtor, cid,
                usage, method, ftype, rcore, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
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
                usage, method, ftype, rcore, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
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
        for icard, card in enumerate(self.cards):
            (eid, caero, aelist, setg, thxi, thyi, dzi, dtori, cidi,
             usagei, methodi, ftypei, rcorei, comment) = card
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
                   coord_id, method, usage, ftype, rcore)
        self.sort()
        self.cards = []

    def _save(self, spline_id, caero_id, aelist_id, set_id,
              thx, thy, dz, dtor,
              coord_id, method, usage, ftype, rcord):
        assert len(self.spline_id) == 0, self.spline_id
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

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.spline_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

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


class GUST(VectorizedBaseCard):
    """
    Defines a stationary vertical gust for use in aeroelastic response
    analysis.

    +------+-----+-------+-----+-----+------+
    |   1  |  2  |   3   |  4  |  5  |  6   |
    +======+=====+=======+=====+=====+======+
    | GUST | SID | DLOAD | WG  | X0  |  V   |
    +------+-----+-------+-----+-----+------+
    | GUST | 133 |   61  | 1.0 | 0.  | 1.+4 |
    +------+-----+-------+-----+-----+------+
    """
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.gust_id = np.array([], dtype='int32')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, sid: int, dload: int, wg: float, x0: float,
            V: Optional[float]=None, comment: str='') -> int:
        """
        Creates a GUST card, which defines a stationary vertical gust
        for use in aeroelastic response analysis.

        Parameters
        ----------
        sid : int
            gust load id
        dload : int
            TLOADx or RLOADx entry that defines the time/frequency
            dependence
        wg : float
            Scale factor (gust velocity/forward velocity) for gust
            velocity
        x0 : float
            Streamwise location in the aerodynamic coordinate system of
            the gust reference point.
        V : float; default=None
            float : velocity of the vehicle (must be the same as the
                    velocity on the AERO card)
            None : ???
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((sid, dload, wg, x0, V, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a GUST card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        dload = integer(card, 2, 'dload')
        wg = double(card, 3, 'wg')
        x0 = double(card, 4, 'x0')
        V = double_or_blank(card, 5, 'V')
        assert len(card) <= 6, f'len(GUST card) = {len(card):d}\ncard={card}'
        #return GUST(sid, dload, wg, x0, V=V, comment=comment)
        self.cards.append((sid, dload, wg, x0, V, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        gust_id = np.zeros(ncards, dtype='int32')
        dload_id = np.zeros(ncards, dtype='int32')

        wg = np.zeros(ncards, dtype='float64')
        x0 = np.zeros(ncards, dtype='float64')
        V = np.zeros(ncards, dtype='float64')

        for icard, card in enumerate(self.cards):
            (gust_idi, dloadi, wgi, x0i, Vi, comment) = card
            gust_id[icard] = gust_idi
            dload_id[icard] = dloadi
            wg[icard] = wgi
            x0[icard] = x0i
            V[icard] = Vi
        self._save(gust_id, dload_id, wg, x0, V)
        #self.sort()
        self.cards = []

    def _save(self, gust_id, dload_id, wg, x0, V):
        if len(self.gust_id):
            gust_id = np.hstack([self.gust_id, gust_id])
            dload_id = np.hstack([self.dload_id, dload_id])
            wg = np.hstack([self.wg, wg])
            x0 = np.hstack([self.x0, x0])
            V = np.hstack([self.V, V])
        self.gust_id = gust_id
        self.dload_id = dload_id
        self.wg = wg
        self.x0 = x0
        self.V = V

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        #geom_check(self,
                   #missing,
                   #coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names))

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.gust_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        gust_ids = array_str(self.gust_id, size=size)
        dload_ids = array_str(self.dload_id, size=size)
        for gust_id, dload, wg, x0, V in zip(gust_ids, dload_ids, self.wg, self.x0, self.V):
            list_fields = ['GUST', gust_id, dload, wg, x0, V]
            bdf_file.write(print_card(list_fields))
        return


class AESTAT(VectorizedBaseCard):
    """
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +--------+------+--------+
    |    1   |   2  |    3   |
    +========+======+========+
    | AESTAT |  ID  | LABEL  |
    +--------+------+--------+
    | AESTAT | 5001 | ANGLEA |
    +--------+------+--------+
    """
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aestat_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, aestat_id: int, label: str, comment: str='') -> int:
        """
        Creates an AESTAT card, which is a variable to be used in a TRIM analysis

        Parameters
        ----------
        aestat_id : int
            unique id
        label : str
            name for the id
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((aestat_id, label, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment='') -> int:
        """
        Adds an AESTAT card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aestat_id = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        assert len(card) <= 3, f'len(AESTAT card) = {len(card):d}\ncard={card}'
        #return AESTAT(aestat_id, label, comment=comment)
        self.cards.append((aestat_id, label, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aestat_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        for icard, card in enumerate(self.cards):
            (aestat_idi, labeli, comment) = card
            aestat_id[icard] = aestat_idi
            label[icard] = labeli
        self._save(aestat_id, label)
        #self.sort()
        self.cards = []

    def _save(self, aestat_id, label):
        if len(self.aestat_id):
            aestat_id = np.hstack([self.aestat_id, aestat_id])
            label = np.hstack([self.label, label])
        self.aestat_id = aestat_id
        self.label = label

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aestat_id) == 0:
            return ''
        #print_card = get_print_card_8_16(size)

        aestat_ids = array_str(self.aestat_id, size=size)
        for aestat_id, label in zip(aestat_ids, self.label):
            list_fields = ['AESTAT', aestat_id, label]
            bdf_file.write(print_card_8(list_fields))
        return


class AEPARM(VectorizedBaseCard):
    """
    Defines a general aerodynamic trim variable degree-of-freedom (aerodynamic
    extra point). The forces associated with this controller will be derived
    from AEDW, AEFORCE and AEPRESS input data.

    +--------+----+--------+-------+
    |    1   | 2  |    3   |   4   |
    +========+====+========+=======+
    | AEPARM | ID | LABEL  | UNITS |
    +--------+----+--------+-------+
    | AEPARM | 5  | THRUST | LBS   |
    +--------+----+--------+-------+
    """
    _id_name = 'aeparm_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aeparm_id = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.units = np.array([], dtype='|U8')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, aeparm_id: int, label: str, units: str,
            comment: str='') -> int:
        """
        Creates an AEPARM card, which defines a new trim variable.

        Parameters
        ----------
        aeparm_id : int
            the unique id
        label : str
            the variable name
        units : str
            unused by Nastran
        comment : str; default=''
            a comment for the card

        """
        self.cards.append((aeparm_id, label, units, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AEPARM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aeparm_id = integer(card, 1, 'aeparm_id')
        label = string(card, 2, 'label')
        units = card.field(3)
        units = '' if units is None else units

        assert len(card) <= 4, f'len(AEPARM card) = {len(card):d}\ncard={card}'
        #return AEPARM(aeparm_id, label, units, comment=comment)
        self.cards.append((aeparm_id, label, units, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aeparm_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        units = np.zeros(ncards, dtype='|U8')

        for icard, card in enumerate(self.cards):
            (aeparm_idi, labeli, unitsi, comment) = card
            aeparm_id[icard] = aeparm_idi
            label[icard] = labeli
            units[icard] = unitsi
        self._save(aeparm_id, label, units)
        self.sort()
        self.cards = []

    def _save(self, aeparm_id, label, units):
        assert len(self.aeparm_id) == 0, self.aeparm_id
        self.aeparm_id = aeparm_id
        self.label = label
        self.units = units

    def __apply_slice__(self, elem: CAERO1, i: np.ndarray) -> None:
        elem.n = len(i)
        elem.aeparm_id = self.aeparm_id[i]
        elem.label = self.label[i]
        elem.units = self.units[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        #geom_check(self,
                   #missing,
                   #coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names))

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aeparm_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        aeparm_ids = array_str(self.aeparm_id, size=size)
        for aeparm_id, label, units in zip(aeparm_ids, self.label, self.units):
            list_fields = ['AEPARM', aeparm_id, label, units]
            bdf_file.write(print_card(list_fields))
        return


class AESURF(VectorizedBaseCard):
    _id_name = 'aesurf_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aesurf_id = np.array([], dtype='int32')

    def add(self, aesurf_id: int, label: str,
            cid1: int, aelist_id1: int,
            cid2: Optional[int]=None, aelist_id2: Optional[int]=None,
            eff: float=1.0, ldw: str='LDW',
            crefc: float=1.0, crefs: float=1.0,
            pllim: float=-np.pi/2., pulim: float=np.pi/2.,
            hmllim=None, hmulim=None, # hinge moment lower/upper limits
            tqllim=None, tqulim=None, # TABLEDi deflection limits vs. dynamic pressure
            comment='') -> int:
        """
        Creates an AESURF card, which defines a control surface

        Parameters
        ----------
        aesid : int
            controller number
        label : str
            controller name
        cid1 / cid2 : int / None
            coordinate system id for primary/secondary control surface
        aelist_id1 / aelist_id2 : int / None
            AELIST id for primary/secondary control surface
        eff : float; default=1.0
            Control surface effectiveness
        ldw : str; default='LDW'
            Linear downwash flag;  ['LDW', 'NODLW']
        crefc : float; default=1.0
            reference chord for the control surface
        crefs : float; default=1.0
            reference area for the control surface
        pllim / pulim : float; default=-pi/2 / pi/2
            Lower/Upper deflection limits for the control surface in radians
        hmllim / hmulim : float; default=None
            Lower/Upper hinge moment limits for the control surface in
            force-length units
        tqllim / tqulim : int; default=None
            Set identification numbers of TABLEDi entries that provide the
            lower/upper deflection limits for the control surface as a
            function of the dynamic pressure
        comment : str; default=''
            a comment for the card

        """
        if cid2 is None:
            cid2 = 0
        if aelist_id2 is None:
            aelist_id2 = 0
        card = (aesurf_id, label, cid1, aelist_id1, cid2, aelist_id2, eff, ldw,
                crefc, crefs, pllim, pulim, hmllim, hmulim,
                tqllim, tqulim, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AESURF card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aesurf_id = integer(card, 1, 'aesid')
        label = string(card, 2, 'label')

        cid1 = integer(card, 3, 'cid1')
        alid1 = integer(card, 4, 'alid1')

        cid2 = integer_or_blank(card, 5, 'cid2', default=0)
        alid2 = integer_or_blank(card, 6, 'alid2', default=0)

        eff = double_or_blank(card, 7, 'eff', default=1.0)
        ldw = string_or_blank(card, 8, 'ldw', default='LDW')
        crefc = double_or_blank(card, 9, 'crefc', default=1.0)
        crefs = double_or_blank(card, 10, 'crefs', default=1.0)

        pllim = double_or_blank(card, 11, 'pllim', default=-np.pi / 2.)
        pulim = double_or_blank(card, 12, 'pulim', default=np.pi / 2.)

        hmllim = double_or_blank(card, 13, 'hmllim')
        hmulim = double_or_blank(card, 14, 'hmulim')
        tqllim = integer_or_blank(card, 15, 'tqllim')
        tqulim = integer_or_blank(card, 16, 'tqulim')
        assert len(card) <= 17, f'len(AESURF card) = {len(card):d}\ncard={card}'
        #return AESURF(aesurf_id, label, cid1, alid1, cid2, alid2, eff, ldw,
                      #crefc, crefs, pllim, pulim, hmllim, hmulim,
                      #tqllim, tqulim, comment=comment)
        card = (aesurf_id, label, cid1, alid1, cid2, alid2, eff, ldw,
                crefc, crefs, pllim, pulim, hmllim, hmulim,
                tqllim, tqulim, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aesurf_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        coord_id = np.zeros((ncards, 2), dtype='int32')
        aelist_id = np.zeros((ncards, 2), dtype='int32')
        eff = np.zeros(ncards, dtype='float64')
        ldw = np.zeros(ncards, dtype='|U3')
        refc = np.zeros(ncards, dtype='float64')
        refs = np.zeros(ncards, dtype='float64')

        pllim = np.zeros(ncards, dtype='float64')
        pulim = np.zeros(ncards, dtype='float64')
        hmllim = np.zeros(ncards, dtype='float64')

        hmulim = np.zeros(ncards, dtype='float64')
        tqllim = np.zeros(ncards, dtype='int32')
        tqulim = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (aesurf_idi, labeli, cid1i, alid1i, cid2i, alid2i, effi, ldwi,
             refci, refsi, pllimi, pulimi, hmllimi, hmulimi,
             tqllimi, tqulimi, commenti) = card
            aesurf_id[icard] = aesurf_idi
            label[icard] = labeli
            coord_id[icard, :] = [cid1i, cid2i]
            aelist_id[icard, :] = [alid1i, alid2i]
            eff[icard] = effi
            ldw[icard] = ldwi
            refc[icard] = refci
            refs[icard] = refsi

            pllim[icard] = pllimi
            pulim[icard] = pulimi

            hmllim[icard] = hmllimi
            hmulim[icard] = hmulimi

            tqllimi = 0 if tqllimi is None else tqllimi
            tqulimi = 0 if tqulimi is None else tqulimi
            tqllim[icard] = tqllimi
            tqulim[icard] = tqulimi

        self._save(aesurf_id, label, coord_id, aelist_id, eff, ldw, refc, refs,
                   pllim, pulim, hmllim, hmulim, tqllim, tqulim)
        self.sort()
        self.cards = []

    def _save(self, aesurf_id, label, coord_id, aelist_id, eff, ldw, refc, refs,
              pllim, pulim, hmllim, hmulim, tqllim, tqulim):
        assert len(self.aesurf_id) == 0, self.aesurf_id
        self.aesurf_id = aesurf_id
        self.label = label
        self.coord_id = coord_id
        self.aelist_id = aelist_id
        self.eff = eff
        self.ldw = ldw
        self.refc = refc
        self.refs = refs

        self.pllim = pllim
        self.pulim = pulim
        self.hmllim = hmllim
        self.hmulim = hmulim
        self.tqllim = tqllim
        self.tqulim = tqulim

        self.n = len(aesurf_id)

    def sort(self) -> None:
        ueid = np.unique(self.aesurf_id)
        if np.array_equal(ueid, self.aesurf_id):
            return
        i = np.argsort(self.aesurf_id)
        self.__apply_slice__(self, i)

    def __apply_slice__(self, elem: AESURF, i: np.ndarray) -> None:
        elem.n = len(i)
        #self.aesurf_id = self.aesurf_id[i]
        #self.label = self.label[i]
        #self.coord_id = self.coord_id[i]
        #self.aelist_id = self.aelist_id[i]
        #self.eff = self.eff[i]
        #self.ldw = self.ldw[i]
        #self.refc = self.refc[i]
        #self.refs = self.refs[i]

        elem.aesurf_id = self.aesurf_id[i]
        elem.label = self.label[i]
        elem.coord_id = self.coord_id[i, :]
        elem.aelist_id = self.aelist_id[i, :]
        elem.eff = self.eff[i]
        elem.ldw = self.ldw[i]
        elem.refc = self.refc[i]
        elem.refs = self.refs[i]

        elem.pllim = self.pllim[i]
        elem.pulim = self.pulim[i]
        elem.hmllim = self.hmllim[i]
        elem.hmulim = self.hmulim[i]
        elem.tqllim = self.tqllim[i]
        elem.tqulim = self.tqulim[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        ucoords = np.unique(self.coord_id)
        uaelist = np.unique(self.aelist_id)
        if uaelist[0] == 0:
            uaelist = uaelist[1:]

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aesurf_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        aesurf_id_ = array_str(self.aesurf_id, size=size)
        cid1_ = array_str(self.coord_id[:, 0], size=size)
        cid2_ = array_str(self.coord_id[:, 1], size=size)

        #elem.label = self.label[i]
        #elem.eff = self.eff[i]
        #elem.ldw = self.ldw[i]
        #elem.refc = self.refc[i]
        #elem.refs = self.refs[i]

        #elem.pllim = self.pllim[i]
        #elem.pulim = self.pulim[i]
        #elem.hmllim = self.hmllim[i]
        #elem.hmulim = self.hmulim[i]
        #elem.tqllim = self.tqllim[i]
        #elem.tqulim = self.tqulim[i]

        aelist_id1_ = array_str(self.aelist_id[:, 0], size=size)
        aelist_id2_ = array_str(self.aelist_id[:, 1], size=size)

        hmllims = array_float_nan(self.hmllim, size=size, is_double=False)
        hmulims = array_float_nan(self.hmulim, size=size, is_double=False)

        tqllims = array_str(self.tqllim, size=size)
        tqulims = array_str(self.tqulim, size=size)

        effs = array_default_float(self.eff, default=0., size=size, is_double=False)
        ldws = array_default_str(self.ldw, default='LDW', size=size)
        for aesurf_id, label, cid1, aelist_id1, \
            cid2, aelist_id2, eff, ldw, crefc, crefs, \
            pllim, pulim, hmllim, hmulim, \
            tqllim, tqulim in zip_longest(
                aesurf_id_, self.label, cid1_, aelist_id1_, cid2_, aelist_id2_,
                effs, ldws, self.refc, self.refs,
                self.pllim, self.pulim, hmllims, hmulims, tqllims, tqulims):
            #eff = set_blank_if_default(eff, 1.0)
            #ldw = set_blank_if_default(ldw, 'LDW')
            crefc = set_blank_if_default(crefc, 1.0)
            crefs = set_blank_if_default(crefs, 1.0)

            pllim = set_blank_if_default(pllim, -np.pi / 2.)
            pulim = set_blank_if_default(pulim, np.pi / 2.)

            list_fields = ['AESURF', aesurf_id, label, cid1, aelist_id1,
                           cid2, aelist_id2, eff, ldw, crefc, crefs,
                           pllim, pulim, hmllim, hmulim, tqllim,
                           tqulim]
            bdf_file.write(print_card(list_fields))
        return


class AESURFS(VectorizedBaseCard):
    """
    Optional specification of the structural nodes associated with an
    aerodynamic control surface that has been defined on an AESURF entry. The
    mass associated with these structural nodes define the control surface
    moment(s) of inertia about the hinge line(s).
    Specifies rigid body motions to be used as trim variables in static
    aeroelasticity.

    +---------+------+-------+---+-------+---+-------+
    |    1    |  2   |   3   | 4 |   5   | 6 |   7   |
    +=========+======+=======+===+=======+===+=======+
    | AESURFS |  ID  | LABEL |   | LIST1 |   | LIST2 |
    +---------+------+-------+---+-------+---+-------+
    | AESURFS | 6001 | ELEV  |   | 6002  |   | 6003  |
    +---------+------+-------+---+-------+---+-------+
    """
    _id_name = 'aesurfs_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.aesurfs_id = np.array([], dtype='int32')

    def add(self, aesurfs_id: int, label: str,
            list1: int, list2: int, comment: str='') -> int:
        """
        Creates an AESURFS card

        Parameters
        ----------
        aesid : int
            the unique id
        label : str
            the AESURF name
        list1 / list2 : int / None
            the list (SET1) of node ids for the primary/secondary
            control surface(s) on the AESURF card
        comment : str; default=''
            a comment for the card

        """
        card = (aesurfs_id, label, list1, list2, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds an AESURFS card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        aesid = integer(card, 1, 'ID')
        label = string(card, 2, 'label')
        list1 = integer(card, 4, 'list1')
        list2 = integer(card, 6, 'list2')
        assert len(card) <= 7, f'len(AESURFS card) = {len(card):d}\ncard={card}'
        #return AESURFS(aesid, label, list1, list2, comment=comment)
        card = (aesid, label, list1, list2, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        aesurfs_id = np.zeros(ncards, dtype='int32')
        label = np.zeros(ncards, dtype='|U8')
        list1_id = np.zeros(ncards, dtype='int32')
        list2_id = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (aesurfs_idi, labeli, list1i, list2i, comment) = card
            aesurfs_id[icard] = aesurfs_idi
            label[icard] = labeli
            list1_id[icard] = list1i
            list2_id[icard] = list2i

        self._save(aesurfs_id, label, list1_id, list2_id)
        self.sort()
        self.cards = []

    def _save(self, aesurfs_id, label, list1_id, list2_id):
        assert len(self.aesurfs_id) == 0, self.aesurfs_id
        self.aesurfs_id = aesurfs_id
        self.label = label
        self.list1_id = list1_id
        self.list2_id = list2_id
        self.n = len(aesurfs_id)

    #def sort(self) -> None:
        #ueid = np.unique(self.aesurf_id)
        #if np.array_equal(ueid, self.aesurf_id):
            #return
        #i = np.argsort(self.aesurf_id)
        #self.__apply_slice__(self, i)

    def __apply_slice__(self, elem: AESURFS, i: np.ndarray) -> None:
        elem.n = len(i)
        #self.aesurf_id = self.aesurf_id[i]
        #self.label = self.label[i]
        #self.coord_id = self.coord_id[i]
        #self.aelist_id = self.aelist_id[i]
        #self.eff = self.eff[i]
        #self.ldw = self.ldw[i]
        #self.refc = self.refc[i]
        #self.refs = self.refs[i]

        elem.aesurfs_id = self.aesurfs_id[i]
        elem.label = self.label[i]
        elem.list1_id = self.list1_id[i]
        elem.list2_id = self.list2_id[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        pass

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.aesurfs_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        aesurfs_id_ = array_str(self.aesurfs_id, size=size)

        list1_id = array_str(self.list1_id, size=size)
        list2_id = array_str(self.list2_id, size=size)

        for aesurfs_id, label, list1, list2 in zip(
                aesurfs_id_, self.label, list1_id, list2_id):

            list_fields = ['AESURFS', aesurfs_id, label, None, list1, None, list2]
            bdf_file.write(print_card(list_fields))
        return


class CSSCHD(VectorizedBaseCard):
    """
    Defines a scheduled control surface deflection as a function of
    Mach number and angle of attack.

    +--------+-----+-------+--------+-------+-------+
    |    1   |  2  |   3   |   4    |   5   |   6   |
    +========+=====+=======+========+=======+=======+
    | CSSCHD | SlD | AESID | LALPHA | LMACH | LSCHD |
    +--------+-----+-------+--------+-------+-------+
    | CSSCHD |  5  |  50   |   12   |   15  |   25  |
    +--------+-----+-------+--------+-------+-------+
    """
    _id_name = 'csschd_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.csschd_id = np.array([], dtype='int32')
        self.aesurf_id = np.array([], dtype='int32')
        self.lalpha = np.array([], dtype='int32')
        self.lmach = np.array([], dtype='int32')
        self.lschd = np.array([], dtype='int32')

    def add(self, sid: int, aesurf_id: int,
            lschd: int, lalpha: int=None, lmach: int=None,  # aefact
            comment: str='') -> int:
        """
        Creates an CSSCHD card, which defines a specified control surface
        deflection as a function of Mach and alpha (used in SOL 144/146).

        Parameters
        ----------
        sid : int
            the unique id
        aesid : int
            the control surface (AESURF) id
        lalpha : int; default=None
            the angle of attack profile (AEFACT) id
        lmach : int; default=None
            the mach profile (AEFACT) id
        lschd : int; default=None
            the control surface deflection profile (AEFACT) id
        comment : str; default=''
            a comment for the card

        """
        #assert lalpha is None or isinstance(lalpha, integer_types), lalpha
        #assert lmach is None or isinstance(lmach, integer_types), lmach
        #assert lschd is None or isinstance(lschd, integer_types), lschd
        card = (sid, aesurf_id, lalpha, lmach, lschd, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a CSSCHD card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        aesurf_id = integer(card, 2, 'aesid')         # AESURF
        lalpha = integer_or_blank(card, 3, 'lAlpha')  # AEFACT
        lmach = integer_or_blank(card, 4, 'lMach')    # AEFACT
        lschd = integer(card, 5, 'lSchd')             # AEFACT
        assert len(card) <= 6, f'len(CSSCHD card) = {len(card):d}\ncard={card}'
        #return CSSCHD(sid, aesurf_id, lalpha, lmach, lschd, comment=comment)
        card = (sid, aesurf_id, lalpha, lmach, lschd, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        csschd_id = np.zeros(ncards, dtype='int32')
        aesurf_id = np.zeros(ncards, dtype='int32')
        lalpha = np.zeros(ncards, dtype='int32')
        lmach = np.zeros(ncards, dtype='int32')
        lschd = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (sid, aesurf_idi, lalphai, lmachi, lschdi, comment) = card
            csschd_id[icard] = sid
            aesurf_id[icard] = aesurf_idi
            lalpha[icard] = lalphai
            lmach[icard] = lmachi
            lschd[icard] = lschdi
        self._save(csschd_id, aesurf_id, lalpha, lmach, lschd)
        self.sort()
        self.cards = []

    def _save(self, csschd_id, aesurf_id, lalpha, lmach, lschd):
        assert len(self.csschd_id) == 0, self.csschd_id
        self.csschd_id = csschd_id
        self.aesurf_id = aesurf_id
        self.lalpha = lalpha
        self.lmach = lmach
        self.lschd = lschd
        self.n = len(csschd_id)

    #def sort(self) -> None:
        #ueid = np.unique(self.csschd_id)
        #if np.array_equal(ueid, self.csschd_id):
            #return
        #i = np.argsort(self.csschd_id)
        #self.__apply_slice__(self, i)

    def __apply_slice__(self, load: CSSCHD, i: np.ndarray) -> None:
        load.n = len(i)
        load.csschd_id = self.csschd_id[i]
        load.aesurf_id = self.aesurf_id[i]
        load.lalpha = self.lalpha[i]
        load.lmach = self.lmach[i]
        load.lschd = self.lschd[i]

    def validate(self):
        #if not(self.lalpha is None or isinstance(self.lalpha, integer_types)):
            #raise TypeError('lalpha=%r must be an int or None' % self.lalpha)
        #if not(self.lmach is None or isinstance(self.lmach, integer_types)):
            #raise TypeError('lmach=%r must be an int or None' % self.lmach)

        ibad = (self.lalpha == 0) and (self.lmach == 0)
        if np.any(ibad):
            csschd = self.csschd_id[ibad]
            card = self.slice_card_by_index(ibad)
            msgi = card.write()
            msg = ('CSSCHD csschd_id=%s; lalpha and lmach are both 0'
                   ' (one must be an integer (AEFACT)\n%s' % (csschd, msgi))
            raise RuntimeError(msg)

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        uaesurf = np.unique(self.aesurf_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.csschd_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        csschd_id_ = array_str(self.csschd_id, size=size)
        aesurf_id_ = array_str(self.aesurf_id, size=size)
        lalpha_ = array_default_int(self.lalpha, default=0, size=size)
        lmach_ = array_default_int(self.lmach, default=0, size=size)
        lschd_ = array_str(self.lschd, size=size)

        for csschd_id, aesurf_id, lalpha, lmach, lschd in zip(
                csschd_id_, aesurf_id_, lalpha_, lmach_, lschd_):

            list_fields = ['CSSCHD', csschd_id, aesurf_id, lalpha, lmach, lschd]
            bdf_file.write(print_card(list_fields))
        return

class TRIM(VectorizedBaseCard):
    """
    Specifies constraints for aeroelastic trim variables.

    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |   1  |   2    |   3  |    4   |    5   |  6  |    7   |  8  |     9    |
    +======+========+======+========+========+=====+========+=====+==========+
    | TRIM |   ID   | MACH |    Q   | LABEL1 | UX1 | LABEL2 | UX2 | IS_RIGID |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    |      | LABEL3 |  UX3 | LABEL4 |   UX4  | ... |        |     |          |
    +------+--------+------+--------+--------+-----+--------+-----+----------+
    """
    _id_name = 'trim_id'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.trim_id = np.array([], dtype='int32')
        self.mach = np.array([], dtype='float64')
        self.q = np.array([], dtype='float64')
        self.aeqr = np.array([], dtype='float64')
        self.nlabel = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U8')
        self.ux = np.array([], dtype='float64')

    def add(self, sid: int, mach: float, q: float,
            labels: list[str], uxs: list[float], aeqr: float=1.0,
            trim_type: int=1, comment: str='') -> int:
        """
        Creates a TRIM/TRIM2 card for a static aero (144) analysis.

        Parameters
        ----------
        sid : int
            the trim id; referenced by the Case Control TRIM field
        mach : float
            the mach number
        q : float
            dynamic pressure
        labels : list[str]
            names of the fixed variables
        uxs : list[float]
            values corresponding to labels
        aeqr : float
            0.0 : rigid trim analysis
            1.0 : elastic trim analysis (default)
        trim_type : int
            1 : creates a TRIM
            2 : creates a TRIM2
        comment : str; default=''
            a comment for the card

        """
        assert len(labels) == len(uxs)
        card = (sid, mach, q, labels, uxs, aeqr, comment)
        self.cards.append(card)
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a TRIM card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        labels = []
        uxs = []

        label = string_or_blank(card, 4, 'label1')
        if label:
            ux = double(card, 5, 'ux1')
            uxs.append(ux)
            labels.append(label)

        label = string_or_blank(card, 6, 'label2')
        if label:
            ux = double(card, 7, 'ux1')
            uxs.append(ux)
            labels.append(label)
        aeqr = double_or_blank(card, 8, 'aeqr', default=1.0)

        i = 9
        n = 3
        while i < len(card):
            label = string(card, i, 'label%i' % n)
            ux = double(card, i + 1, 'ux%i' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
            n += 1
        #return TRIM(sid, mach, q, labels, uxs, aeqr, comment=comment)
        self.cards.append((sid, mach, q, labels, uxs, aeqr, comment))
        self.n += 1
        return self.n - 1

    def add_card2(self, card: BDFCard, comment: str='') -> int:
        """
        Adds a TRIM2 card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card

        """
        sid = integer(card, 1, 'sid')
        mach = double(card, 2, 'mach')
        q = double(card, 3, 'q')
        aeqr = double_or_blank(card, 8, 'aeqr', default=1.0)

        i = 9
        n = 3
        labels = []
        uxs = []
        while i < len(card):
            label = string(card, i, 'label%d' % n)
            ux = double(card, i + 1, 'ux%d' % n)
            labels.append(label)
            uxs.append(ux)
            i += 2
        #return TRIM2(sid, mach, q, labels, uxs, aeqr, comment=comment)
        self.cards.append((2, sid, mach, q, labels, uxs, aeqr, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        trim_id = np.zeros(ncards, dtype='int32')
        mach = np.zeros(ncards, dtype='float64')
        q = np.zeros(ncards, dtype='float64')
        aeqr = np.zeros(ncards, dtype='float64')
        nlabel = np.zeros(ncards, dtype='int32')
        label = []
        ux = []

        for icard, card in enumerate(self.cards):
            (sid, machi, qi, labelsi, uxsi, aeqri, comment) = card
            trim_id[icard] = sid
            mach[icard] = machi
            q[icard] = qi
            aeqr[icard] = aeqri
            nlabel[icard] = len(labelsi)
            label.extend(labelsi)
            ux.extend(uxsi)
        label = np.array(label, dtype='|U8')
        ux = np.array(ux, dtype='float64')
        self._save(trim_id, mach, q, aeqr, nlabel, label, ux)
        self.sort()
        self.cards = []

    def _save(self, trim_id, mach, q, aeqr, nlabel, label, ux):
        assert len(self.trim_id) == 0, self.trim_id
        self.trim_id = trim_id
        self.mach = mach
        self.q = q
        self.aeqr = aeqr
        self.nlabel = nlabel
        self.label = label
        self.ux = ux
        self.n = len(trim_id)

    def validate(self):
        assert self.mach.min() >= 0.0, 'mach = %s' % self.mach
        assert np.all(self.mach != 1.0), 'mach = %s' % self.mach
        assert self.q.min() > 0.0, 'q=%s' % self.q

        for trim_id, (ilabel0, ilabel1) in zip(self.trim_id, self.ilabel):
            labels = self.label[ilabel0:ilabel1]
            if len(set(labels)) != len(labels):
                msg = 'TRIM id=%d; not all labels are unique; labels=%s' % (trim_id, str(labels))
                raise RuntimeError(msg)
        #if len(self.labels) != len(self.uxs):
            #msg = 'nlabels=%s != nux=%s; labels=%s uxs=%s' % (
                #len(self.labels), len(self.uxs), str(self.labels), str(self.uxs))
            #raise RuntimeError(msg)

    @property
    def ilabel(self) -> np.ndarray:
        return make_idim(self.n, self.nlabel)

    def __apply_slice__(self, load: TRIM, i: np.ndarray) -> None:
        load.n = len(i)
        load.trim_id = self.trim_id[i]
        load.mach = self.mach[i]
        load.q = self.q[i]
        load.lmach = self.lmach[i]

        ilabel = self.ilabel
        load.label = hslice_by_idim(i, ilabel, self.label)
        load.ux = hslice_by_idim(i, ilabel, self.ux)
        load.nlabel = self.nlabel[i]

    def geom_check(self, missing: dict[str, np.ndarray]):
        model = self.model
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        #coords = self.model.coord.coord_id
        uaesurf = np.unique(self.trim_id)

        #set1_ids = np.unique(set1_ids)
        geom_check(
            self,
            missing,
            #coord=(model.coord.coord_id, ucoords),
            #aelist=(model.aelist.aelist_id, aelist_ids),
            #caero=(model.caero1.caero_id, caero_ids),
        )

    def verify_trim(self, suport1_id: int) -> None:
        """
        Magic function that makes TRIM cards not frustrating.

        .. warning ::  This probably gets AELINKs/AEPARMs/AESURFSs wrong.

        **The TRIM equality**
        ndelta = (naestat + naesurf + naeparm)
               - (ntrim + ntrim_aesurf? + naelink + nsuport_dofs + nsuport1_dofs)
        ndelta = 0
        ntrim_aesurf is not included, but it might exist...

        **Steps to a TRIM analysis**
        1.  Define the number of independent control surfaces (naesurf)
            Make an AESURF for each.  Dual link the AESURFs if you can
            to avoid needing an AELINK (e.g., +roll is left aileron down,
            right aileron up).
            Horizontal Tail : name it DPITCH
            Vertical Tail   : name it DYAW
            Aileron         : name it DROLL
        2.  Create AELINKs if necessary.
        3.  Add the AESTAT variables.  Include one for each DOF the
            aircraft can move in the frame of the model
            (e.g., half/full model).
                Half model (2.5g pullup, abrupt pitch):
                  - 2d pitch/plunge, 1 control : URDD3, URDD5, PITCH, ANGLEA
                Full model (2.5g pullup, abrupt pitch):
                  - 3d pitch/plunge, 3 control : URDD3, URDD5, PITCH, ANGLEA, YAW (???)
        4.  Add the TRIM card to lock the variables that could theoretically move
            in the plane of the analysis that are known.
                Half model:
                   2.5g pullup   : lock URDD3=2.5, URDD5=0, PITCH=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH
                   abrupt pitch  : lock URDD3=1.0, URDD5=0, ANGLEA=5
                                   solve for PITCH, DPITCH
                                   use DPITCH
                Full model:
                   2.5g pullup   : lock URDD3=2.5, URDD4=0, URDD5=0,  PITCH=0, YAW=0,
                                   lock SIDES=0,  ROLL=0
                                   solve for ANGLEA, DPITCH
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong
                   30 degree yaw : lock URDD3=1.0, URDD4=0, ANGLEA=5, PITCH=0, YAW=30,
                                   lock DPITCH=0, ROLL=0
                                   solve for SIDES, URDD5
                                   use DPITCH, DYAW, DROLL
                                   TODO: probably wrong

        5.  Note that we could have simplified our full model AESTAT/TRIM
            cards (they can be the same as for a half model), but we'd
            like to be able to do multiple load cases in the same deck.

        6.  Add some SUPORT/SUPORT1 DOFs to ignore non-relevant motion in
            certain DOFs (e.g., z-motion).  Add enough to satisfy the TRIM
            equality.

        **Doesn't Consider**
         - AELINK
         - AEPARM
         - AESURFS

        +------------------------------------------------+
        |                 Default AESTATs                |
        +--------+---------+-----------------------------+
        | ANGLEA | ur (R2) | Angle of Attack             |
        | YAW    | ur (R3) | Yaw Rate                    |
        | SIDES  | ur (R3) | Angle of Sideslip           |
        +--------+---------+-----------------------------+
        | ROLL   | ůr (R1) | Roll Rate                   |
        | PITCH  | ůr (R2) | Pitch Rate                  |
        +--------+---------+-----------------------------+
        | URDD1  | ür (T1) | Longitudinal (See Remark 3) |
        | URDD2  | ür (T2) | Lateral                     |
        | URDD3  | ür (T3) | Vertical                    |
        | URDD4  | ür (R1) | Roll                        |
        | URDD5  | ür (R2) | Pitch                       |
        | URDD6  | ür (R3) | Yaw                         |
        +--------+---------+-----------------------------+
        """
        #if not xref:
            #return
        #suport, suport1, aestats, aeparms, aelinks, aesurf, xref=True
        #suport = []
        #if 'SUPORT1' in subcase:
            #suport_id = subcase.get_int_parameter('SUPORT1')[0]
            #suport1 = fem.suport1[suport_id]

        suport = self.model.suport

        suport_ids = np.unique(suport.suport_id)
        nsuport_dofs = 0
        nsuport1_dofs = 0
        suport_dofs = set()

        #for (inode0, inode1) in suport.inode:
            #nodes = suport.node[inode0:inode1]
            #components = suport.component[inode0:inode1].astype('|U8')
            #for nid, cs in zip(nodes, components):
                #for ci in str(cs):
                    ##print('  SUPORT: nid=%r C=%r' % (nid, ci))
                    #dof = (nid, ci)
                    #if dof in suport_dofs:
                        #msg = 'Duplicate DOF\n  dof=%s suport_dofs=%s' % (
                            #str(dof), str(suport_dofs))
                        #raise RuntimeError(msg)
                    #suport_dofs.add(dof)
                    #nsuport_dofs += 1

        suport_dof_msg2 = ''
        if 0 in suport_ids:
            suport_id = 0
            suport0 = suport.slice_card_by_id(suport_id) # , assume_sorted=True, sort_ids=False
            suport_dof_msg = ''
            for nid, component in zip(suport0.node_id, suport0.component):
                for componenti in str(component):
                    dof = (nid, componenti)
                    suport_dof_msg += '    (%s, %s)\n' % (nid, componenti)
                    if dof in suport_dofs:
                        msg = 'dof=%s suport_dofs=%s' % (str(dof), str(suport_dofs))
                        raise RuntimeError(msg)
                    suport_dofs.add(dof)
                    nsuport_dofs += 1
            suport_dof_msg2 = '\nsuport_dofs (nid, comp):\n%s\n' % suport_dof_msg.rstrip(',')

        if suport1_id > 0:
            suport1 = suport.slice_card_by_id(suport1_id) # , assume_sorted=True, sort_ids=False
            suport1_dof_msg = ''
            for nid, component in zip(suport1.node_id, suport1.component):
                for componenti in str(component):
                    dof = (nid, componenti)
                    suport1_dof_msg += '    (%s, %s)\n' % (nid, componenti)
                    if dof in suport1_dofs:
                        msg = 'dof=%s suport1_dofs=%s' % (str(dof), str(suport_dofs))
                        raise RuntimeError(msg)
                    suport_dofs.add(dof)
                    nsuport1_dofs += 1
            suport_dof_msg2 = '\nsuport_dofs (nid, comp):\n%s\n' % suport1_dof_msg.rstrip(',')

        aesurf_names = self.model.aesurf.label.tolist()
        aestat_labels = self.model.aestat.label.tolist()
        aeparm_labels = self.model.aeparm.label.tolist()
        print(self.get_stats())
        #aesurf_names = [aesurfi.label for aesurfi in aesurf.values()]
        #aestat_labels = [aestat.label for aestat in aestats.values()]
        #aeparm_labels = [aeparm.label for aeparm in aeparms.values()]
        naestat = len(aestat_labels)
        ntrim = len(self.label)
        trim_aesurf_common = list(set(self.label).intersection(set(aesurf_names)))
        trim_aesurf_common.sort()
        ntrim_aesurfs = len(trim_aesurf_common)
        naesurf = len(aesurf_names)
        naeparm = len(aeparm_labels)

        aelinksi = []
        trim_id = self.trim_id[0]
        (ilabel0, ilabel1) = self.ilabel[0]
        labels = self.label[ilabel0:ilabel1]

        i0 = np.where((self.model.aelink.aelink_id == 0) |
                       self.model.aelink.aelink_id == trim_id)[0]
        if len(i0):
            aelinksi = self.model.aelink.slice_card_by_index(i0)
        #if 'ALWAYS' in aelinks:
            #aelinksi += [aelink.label for aelink in aelinks['ALWAYS']]

        naelink = len(aelinksi)


        ntrim_aesurf = 0
        allowed_labels = aestat_labels + aesurf_names + aeparm_labels
        msg = ''
        for label in labels:
            if label not in allowed_labels:
                msg += 'TRIM label=%r is not defined\n' % label

            if label in aesurf_names:
                #print('AESTAT/AESURF label = %r' % label)
                ntrim_aesurf += 1
        if msg:
            msg += '\n aestat_labels=%s\n aeparm_labels=%s\n aesurf_names=%s\n%s' % (
                aestat_labels, aeparm_labels, aesurf_names, str(self))
            raise RuntimeError(msg)

        # TODO: this doesn't work for multiple subcases
        #ntotal_suport_dofs = nsuport_dofs, nsuport1_dofs
        #ndelta = ntrim - nsuport_dofs - nsuport1_dofs - naesurf
        #if ndelta != 0:
            #msg = 'ntrim - nsuport_dofs - nsuport1_dofs - naesurf = ndelta = %s; ndelta != 0\n' % ndelta
            #msg += 'ntrim=%s nsuport_dofs=%s nsuport1_dofs=%s naesurfs=%s' % (
                #ntrim, nsuport_dofs, nsuport1_dofs, naesurf)
            #raise RuntimeError(msg)

        #ndelta = (naestat + naesurf + naeparm + ntrim_aesurf) - (ntrim + naelink + nsuport_dofs + nsuport1_dofs)
        #if ndelta != 0:
            #msg = (
                #'(naestat + naesurf + naeparm + ntrim_aesurf) - '
                #'(ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                #'naestat=%s naesurf=%s naeparm=%s ntrim_aesurfs=%s\n'
                #'ntrim=%s naelink=%s nsuport_dofs=%s nsuport1_dofs=%s' % (
                    #ndelta,
                    #naestat, naesurf, naeparms, ntrim_aesurf,
                    #ntrim, naelink, nsuport_dofs, nsuport1_dofs))

        nplus = (naestat + naesurf + naeparm)
        nminus = ntrim + naelink + nsuport_dofs + nsuport1_dofs

        ndelta = nplus - nminus + 0*2*ntrim_aesurfs
        if ndelta != 0:
            #msg = (
                #'(naestat + naesurf + naeparm) - (ntrim + ntrim_aesurf? + naelink + '
                #'nsuport_dofs + nsuport1_dofs) = ndelta = %s; ndelta != 0\n'
                #'naestat=%s naesurf=%s naeparm=%s ntrim=%s ntrim_aesurf=%s '
                #'naelink=%s nsuport_dofs=%s nsuport1_dofs=%s\n' % (
                    #ndelta,
                    #naestat, naesurf, naeparm, ntrim, ntrim_aesurf,
                    #naelink, nsuport_dofs, nsuport1_dofs)
            #)
            msg = (
                'Invalid trim state (ndelta != 0):\n'
                f'   (naestat + naesurf + naeparm + 0*2*ntrim_aesurf?) = ({naestat} + {naesurf} + {naeparm} + 0*2*{ntrim_aesurf}) = {nplus}\n'
                f' - (ntrim + naelink + nsuport_dofs + nsuport1_dofs) = ({ntrim} + {naelink} + {nsuport_dofs} + {nsuport1_dofs}) = {nminus}\n'
                '===================================================================\n'
                f'  ndelta = {ndelta}\n\n'
                'Summary\n'
                '-------\n'
                f'  +naestat = {naestat}; {aestat_labels}\n'
                f'  +naesurf = {naesurf}; {aesurf_names}\n'
                f'  +naeparm = {naeparm}; {aeparm_labels}\n'
                f'  +0*2*ntrim_aesurf? = {2*ntrim_aesurf} -> 0; {trim_aesurf_common}\n'
                f'  -ntrim = {ntrim}; {labels}\n'
                f'  -naelink = {naelink}; {aelinksi}\n'
                f'  -nsuport_dofs = {nsuport_dofs}\n'
                f'  -nsuport1_dofs = {nsuport1_dofs}\n'
                f'{suport_dof_msg2}\n\n'
            )
            msg += str(self)
            raise RuntimeError(msg)

    def write_file(self, bdf_file: TextIOLike, size: int=8,
                   is_double: bool=False,
                   write_card_header: bool=False) -> None:
        if len(self.trim_id) == 0:
            return ''
        print_card = get_print_card_8_16(size)

        trim_ids = array_str(self.trim_id, size=size)

        for trim_id, mach, q, aeqr, (ilabel0, ilabel1) in zip(
                trim_ids, self.mach, self.q, self.aeqr, self.ilabel):

            labels = self.label[ilabel0:ilabel1]
            uxs = self.ux[ilabel0:ilabel1]
            list_fields = ['TRIM', trim_id, mach, q]
            nlabels = len(labels)
            assert nlabels > 0, labels
            for (i, label, ux) in zip(count(), labels, uxs):
                list_fields += [label, ux]
                if i == 1:
                    list_fields += [aeqr]
            if nlabels == 1:
                list_fields += [None, None, aeqr]
            bdf_file.write(print_card(list_fields))
        return
