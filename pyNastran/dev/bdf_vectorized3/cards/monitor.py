from __future__ import annotations
#from abc import abstractmethod
#from itertools import count
from typing import Optional, TYPE_CHECKING
import numpy as np

from pyNastran.utils.numpy_utils import integer_types # , cast_ints
#from pyNastran.bdf.field_writer_8 import print_card_8, set_blank_if_default
#from pyNastran.bdf.field_writer_16 import print_card_16 # , print_scientific_16, print_field_16
#from pyNastran.bdf.field_writer_double import print_scientific_double
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, # integer_or_string,
    double_or_blank,
    string, string_or_blank, parse_components,
    components_or_blank,
)
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
#from pyNastran.bdf.cards.elements.bars import set_blank_if_default
#from pyNastran.bdf.cards.base_card import expand_thru

from pyNastran.dev.bdf_vectorized3.cards.base_card import (
    VectorizedBaseCard, parse_check) # , make_idim, hslice_by_idim
from pyNastran.dev.bdf_vectorized3.cards.write_utils import array_str, array_default_int, array_float, get_print_card_size
from pyNastran.dev.bdf_vectorized3.bdf_interface.geom_check import geom_check
#from pyNastran.femutils.utils import hstack_lists


if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
    #from pyNastran.nptyping_interface import NDArray3float
    #from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike

ELEMENT_TABLES = {'FORCE', 'STRESS', 'STRAIN'}
ELEMENT_TYPES = {
    #'CELAS1', 'CELAS2', 'CELAS3',
    #'CBUSH', 'CBUSH1D', 'CGAP',
    #'CONROD', 'CROD', 'CTUBE',
    #'CBAR', 'CBEAM', 'CBEND',
    #'CQUAD4', 'CTRIA3', 'CSHEAR',
    'BEAM', 'QUAD4', # 'CQUAD144',
    'CQUAD4',
    'CQUADR', 'CQUADX4', 'QUAD144',
    #'CHEXA', 'CPYRAM', 'CPENTA', 'CTETRA',
    #'CPLSTN3', 'CPLSTN4', 'CPLSTN6', 'CPLSTN8',
    #'CPLSTS3', 'CPLSTS4', 'CPLSTS6', 'CPLSTS8',
}
class MONPNT1(VectorizedBaseCard):
    """
    +---------+---------+------+-----+-----+-------+------+----+----+
    |    1    |    2    |  3   |  4  |  5  |   6   |   7  | 8  | 9  |
    +=========+=========+======+=====+=====+=======+======+====+====+
    | MONPNT1 |  NAME   |                   LABEL                   |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |  AXES   | COMP | CP  |  X  |   Y   |   Z  | CD |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    | MONPNT1 | WING155 |    Wing Integrated Load to Butline 155    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    |         |    34   | WING |     | 0.0 | 155.0 | 15.0 |    |    |
    +---------+---------+------+-----+-----+-------+------+----+----+
    """
    _id_name = 'name'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U72')
        #self.table = np.array([], dtype='|U8')
        self.axes = np.array([], dtype='|U8')
        self.comp = np.array([], dtype='|U8')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self.cp = np.array([], dtype='int32')
        self.cd = np.array([], dtype='int32')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, name: str, label: str, axes: str, aecomp_name: str,
            xyz: list[float],
            cp: int=0,
            cd: Optional[int]=None, comment: str='') -> int:
        """
        Creates a MONPNT1 card

        Parameters
        ----------
        name : str
            Character string of up to 8 characters identifying the
            monitor point
        label : str
            A string comprising no more than 56 characters
            that identifies and labels the monitor point.
        axes : str
            components {1,2,3,4,5,6}
        aecomp_name : str
            name of the AECOMP/AECOMPL entry
        xyz : list[float, float, float]; default=None
            The coordinates in the CP coordinate system about which the
            loads are to be monitored.
            None : [0., 0., 0.]
        cp : int, CORDx; default=0
           int : coordinate system
        cd : int; default=None -> cp
            the coordinate system for load outputs
        comment : str; default=''
            a comment for the card

        Notes
        -----
        CD - MSC specific field

        """
        if cd is None:
            cd = -1
        self.cards.append((name, label, axes, aecomp_name, cp, xyz, cd, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        name = string(card, 1, 'name')
        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        comp = string(card, 10, 'comp')
        cp = integer_or_blank(card, 11, 'cp', default=0)
        xyz = [
            double_or_blank(card, 12, 'x', default=0.0),
            double_or_blank(card, 13, 'y', default=0.0),
            double_or_blank(card, 14, 'z', default=0.0),
        ]
        cd = integer_or_blank(card, 15, 'cd', default=cp)
        self.cards.append((name, label, axes, comp, cp, xyz, cd, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        label = np.zeros(ncards, dtype='|U72')
        comp = np.zeros(ncards, dtype='|U8')
        axes = np.zeros(ncards, dtype='int32')
        cp = np.zeros(ncards, dtype='int32')
        xyz = np.zeros((ncards, 3), dtype='float64')
        cd = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (namei, labeli, axesi, compi, cpi, xyzi, cdi, comment) = card
            name[icard] = namei
            label[icard] = labeli
            comp[icard] = compi
            axes[icard] = axesi
            cp[icard] = cpi
            xyz[icard, :] = xyzi
            cd[icard] = cdi
        self._save(name, label, comp, axes, cp, xyz, cd)
        self.sort()
        self.cards = []

    def _save(self, name, label, comp, axes, cp, xyz, cd):
        assert len(self.name) == 0, self.name
        self.name = name
        self.label = label
        self.comp = comp
        self.axes = axes
        self.cp = cp
        self.xyz = xyz
        self.cd = cd

    def __apply_slice__(self, monitor: MONPNT1, i: np.ndarray) -> None:  # ignore[override]
        monitor.n = len(i)
        monitor.name = self.name[i]
        monitor.label = self.label[i]
        monitor.comp = self.comp[i]
        monitor.axes = self.axes[i]
        monitor.cp = self.cp[i]
        monitor.xyz = self.xyz[i, :]
        monitor.cd = self.cd[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        # comp
        used_dict['coord_id'].append(self.cp)
        cd = self.cd[self.cd >= 0]
        used_dict['coord_id'].append(cd)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        ucoords = np.unique(np.hstack([self.cp, self.cd]))
        geom_check(self,
                   missing,
                   coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names),
                   )

    @property
    def max_id(self) -> int:
        return max(self.cp.max(), self.cd.max())

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
              size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> None:
        size = 8
        print_card, size = get_print_card_size(size, self.max_id)
        assert size == 8, size

        #cp = self.Cp()
        #cd = self.Cd()

        ## Default = the coordinate system specified by the CP field
        #if cd == cp:
            #cd = ''

        cps = array_default_int(self.cp, size=size)
        cds = array_str(self.cd, size=size)
        is_default_cd = (self.cp == self.cd)
        cds[is_default_cd] = ''
        #cds = array_default_int(self.cd, default=0, size=size)
        xyzs = array_float(self.xyz, size=size, is_double=False).tolist()
        for name, label, comp, axes, cp, xyz, cd in zip(self.name, self.label, self.comp, self.axes,
                                                        cps, xyzs, cds):
            x, y, z = xyz
            msg = 'MONPNT1 %-8s%s\n' % (name, label)
            msg += '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s\n' % (
                axes, comp, cp, x, y, z,
                cd)
            bdf_file.write(msg)
        return


class MONPNT2(VectorizedBaseCard):
    """NX/MSC Nastran card"""
    _id_name = 'name'
    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='int32')
        self.label = np.array([], dtype='|U72')
        self.table = np.array([], dtype='|U8')
        self.element_type = np.array([], dtype='|U8')
        self.nddl_item = np.array([], dtype='int32')
        self.element_id = np.array([], dtype='int32')

    def add(self, name: str, label: str, table: str, element_type: str,
            nddl_item: int, eid: int, comment: str='') -> int:
        """Creates a MONPNT2 card

         Parameters
         ----------
         nddl_item : str (nx)
             ???
         """
        self.cards.append((name, label, table, element_type, nddl_item, eid, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        table = string(card, 9, 'table')
        element_type = string(card, 10, 'type')
        nddl_item = string(card, 11, 'comp/nddl_item')  #  nx = string
        #nddl_item = integer_or_string(card, 11, 'comp/nddl_item')
        #nddl_item = integer_or_blank(card, 11, 'nddl_item')
        eid = integer_or_blank(card, 12, 'eid')
        #return MONPNT2(name, label, table, Type, nddl_item, eid, comment=comment)
        self.cards.append((name, label, table, element_type, nddl_item, eid, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        """
        Table   Item
        ------  ------
        STRESS  CQUAD4
        """
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        label = np.zeros(ncards, dtype='|U72')
        table = np.zeros(ncards, dtype='|U8')
        element_type = np.zeros(ncards, dtype='|U8')
        nddl_item = np.zeros(ncards, dtype='|U8')
        element_id = np.zeros(ncards, dtype='int32')

        #table_int = np.zeros(ncards, dtype='int32')
        #nddl_item_int = np.zeros(ncards, dtype='int32')

        for icard, card in enumerate(self.cards):
            (namei, labeli, tablei, element_typei, nddl_itemi, eid, comment) = card
            name[icard] = namei
            label[icard] = labeli
            table[icard] = tablei
            element_type[icard] = element_typei
            assert tablei in ELEMENT_TABLES, f'table={tablei}; card={card}'
            assert element_typei in ELEMENT_TYPES, f'element_type={element_typei}; card={card}'
            nddl_item[icard] = nddl_itemi
            element_id[icard] = eid
        self._save(name, label, table, element_type, nddl_item, element_id)
        self.sort()
        self.cards = []

    def _save(self, name, label, table, element_type, nddl_item, element_id):
        assert len(self.name) == 0, self.name
        self.name = name
        self.label = label
        self.table = table
        self.element_type = element_type
        self.nddl_item = nddl_item
        self.element_id = element_id

    def __apply_slice__(self, monitor: MONPNT2, i: np.ndarray) -> None:  # ignore[override]
        monitor.n = len(i)
        monitor.name = self.name[i]
        monitor.label = self.label[i]
        monitor.table = self.table[i]
        monitor.element_type = self.element_type[i]
        monitor.nddl_item = self.nddl_item[i]
        monitor.element_id = self.element_id[i]

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['element_id'].append(self.element_id)

    #def geom_check(self, missing: dict[str, np.ndarray]):
        ##mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          ##msg=f'no materials for {self.type}')
        ##mids.sort()
        #coords = self.model.coord.coord_id
        ##all_aecomp_names = self.model.aecomp.name
        ##aecomp_names = np.unique(self.comp)
        #ucoords = np.unique(np.hstack([self.cp, self.cd]))
        #geom_check(self,
                   #missing,
                   #coord=(coords, ucoords),
                   ##aecomp=(all_aecomp_names, aecomp_names),
                   #)

    @property
    def max_id(self) -> int:
        return self.element_id.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
              size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> None:
        size = 8
        print_card, size = get_print_card_size(size, self.max_id)
        assert size == 8, size

        for name, label, table, element_type, nddl_item, eid in zip(self.name, self.label, self.table, self.element_type,
                                                            self.nddl_item, self.element_id):
            msg = 'MONPNT2 %-8s%s\n' % (name, label)
            msg += ('        %-8s%-8s%-8s%-8s\n' % (
                table, element_type, nddl_item, eid
            ))
            #card = self.repr_fields()
            #return self.comment + msg.rstrip() + '\n'
            bdf_file.write(msg)
        return


class MONPNT3(VectorizedBaseCard):
    """MSC Nastran specific card"""
    _id_name = 'name'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='|U8')
        self.label = np.array([], dtype='|U72')
        self.axes = np.array([], dtype='int32')
        self.grid_set = np.array([], dtype='int32')
        self.elem_set = np.array([], dtype='int32')
        self.cp = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self.cd = np.array([], dtype='int32')
        self.xflag = np.array([], dtype='|U1')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, name: str, label: str, axes: str,
            grid_set: int, xyz: list[float],
            elem_set: int=0,
            cp: int=0, cd: Optional[int]=None,
            xflag='', comment: str='') -> int:
        """Creates a MONPNT3 card"""
        if cd is None:
            # ew can do this because this is MSC-specific
            cd = cp
        elem_set = 0 if elem_set is None else elem_set
        assert isinstance(elem_set, integer_types), elem_set
        self.cards.append((name, label, axes, grid_set, elem_set, xyz,
                           cp, cd, xflag, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        name = string(card, 1, 'name')

        label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        label = ''.join(label_fields).strip()
        assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        grid_set = integer(card, 10, 'grid_set')
        elem_set = integer_or_blank(card, 11, 'elem_set', default=0)
        cp = integer_or_blank(card, 12, 'cp', default=0)
        xyz = [
            double_or_blank(card, 13, 'x', default=0.0),
            double_or_blank(card, 14, 'y', default=0.0),
            double_or_blank(card, 15, 'z', default=0.0),
        ]
        xflag = string_or_blank(card, 16, 'xflag', default='')
        cd = integer_or_blank(card, 17, 'cd', default=cp)
        #return MONPNT3(name, label, axes, grid_set, elem_set, xyz,
                       #cp=cp, cd=cd, xflag=xflag, comment=comment)
        self.cards.append((name, label, axes, grid_set, elem_set, xyz,
                           cp, cd, xflag, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        label = np.zeros(ncards, dtype='|U72')
        axes = np.zeros(ncards, dtype='int32')
        grid_set = np.zeros(ncards, dtype='int32')
        elem_set = np.zeros(ncards, dtype='int32')
        cp = np.zeros(ncards, dtype='int32')
        xyz = np.zeros((ncards, 3), dtype='float64')
        cd = np.zeros(ncards, dtype='int32')
        xflag = np.zeros(ncards, dtype='|U1')
        for icard, card in enumerate(self.cards):
            (namei, labeli, axesi, grid_seti, elem_seti, xyzi,
             cpi, cdi, xflagi, comment) = card
            name[icard] = namei
            label[icard] = labeli
            axes[icard] = axesi
            grid_set[icard] = grid_seti
            elem_set[icard] = elem_seti
            cp[icard] = cpi
            xyz[icard, :] = xyzi
            cd[icard] = cdi
            xflag[icard] = xflagi
        self._save(name, label, axes, grid_set, elem_set, cp, xyz, cd, xflag)
        #self.sort()
        self.cards = []

    def _save(self, name, label, axes, grid_set, elem_set, cp, xyz, cd, xflag):
        assert len(self.name) == 0, self.name
        self.name = name
        self.label = label
        self.axes = axes
        self.grid_set = grid_set
        self.elem_set = elem_set
        self.cp = cp
        self.xyz = xyz
        self.cd = cd
        self.xflag = xflag

    def __apply_slice__(self, monitor: MONPNT3, i: np.ndarray) -> None:  # ignore[override]
        monitor.name = self.name[i]
        monitor.label = self.label[i]
        monitor.axes = self.axes[i]
        monitor.grid_set = self.grid_set[i]
        monitor.elem_set = self.elem_set[i]
        monitor.cp = self.cp[i]
        monitor.xyz = self.xyz[i, :]
        monitor.cd = self.cd[i]
        monitor.xflag = self.xflag[i]
        monitor.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['coord_id'].append(self.cp)
        used_dict['coord_id'].append(self.cd)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        ucoords = np.unique(np.hstack([self.cp, self.cd]))
        geom_check(self,
                   missing,
                   coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names),
                   )

    @property
    def max_id(self) -> int:
        return self.grid_set.max()

    @parse_check
    def write_file(self, bdf_file: TextIOLike,
              size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> None:
        size = 8
        print_card, size = get_print_card_size(size, self.max_id)
        assert size == 8, size


        #cp = self.Cp()
        #cd = self.Cd()

        ## Default = the coordinate system specified by the CP field
        #if cd == cp:
            #cd = ''

        cps = array_default_int(self.cp, size=size)
        cds = array_str(self.cd, size=size)
        grid_sets = array_str(self.grid_set, size=size)
        elem_sets = array_str(self.elem_set, size=size)
        is_default_cd = (self.cp == self.cd)
        cds[is_default_cd] = ''
        #cds = array_default_int(self.cd, default=0, size=size)
        xyzs = array_float(self.xyz, size=size, is_double=False).tolist()
        for name, label, axes, grid_set, elem_set, \
            cp, xyz, xflag, cd in zip(self.name, self.label, self.axes,
                                      grid_sets, elem_sets,
                                      cps, xyzs, self.xflag, cds):
            x, y, z = xyz
            msg = (
                'MONPNT3 %-8s%s\n'
                '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n'
                '        %-8s' % (
                    name, label,
                    axes, grid_set, elem_set, cp,
                    x, y, z, xflag,
                    cd))

            #msg = 'MONPNT3 %-8s%s\n' % (self.name, self.label)
            #msg += ('        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n'
                    #'         %-8s' % (
                        #self.axes, self.grid_set, self.elem_set, cp,
                        #print_float_8(x), print_float_8(y), print_float_8(z),
                        #xflag, cd
                    #))
            bdf_file.write(msg.rstrip('\n ') + '\n')
        return


class MONDSP1(VectorizedBaseCard):
    _id_name = 'name'

    @VectorizedBaseCard.clear_check
    def clear(self) -> None:
        self.name = np.array([], dtype='|U8')
        self.label = np.array([], dtype='|U72')
        self.axes = np.array([], dtype='int32')
        self.component = np.array([], dtype='int32')
        self.ind_dof = np.array([], dtype='int32')
        self.cp = np.array([], dtype='int32')
        self.xyz = np.zeros((0, 3), dtype='float64')
        self.cd = np.array([], dtype='int32')

    #def __len__(self) -> int:
        #return len(self.name)

    def add(self, name: str, label: str, axes: str,
            component: int, xyz: list[float],
            cp: int=0, cd: Optional[int]=None,
            ind_dof: str='123', comment: str='') -> int:
        """Creates a MONDSP1 card"""
        cd = cp if cd is None else cd
        self.cards.append((name, label, axes, component, xyz, cp,
                           cd, ind_dof, comment))
        self.n += 1
        return self.n - 1

    def add_card(self, card: BDFCard, comment: str='') -> int:
        row0 = card[0]
        row1 = card[1]
        assert len(card) == 2, card
        assert len(row0) > 8, row0
        assert ',' not in row1, row1
        if '\t' in row1:
            card_fields = row1.split('\t')
            row1 = row1.expandtabs(tabsize=8)

        name = row0[8:16]
        label = row0[16:72]
        card_fields = [
            'MONDSP1', name,
            label[:8], label[8:16], label[16:24], label[24:32], label[32:40], label[40:48], label[48:56],
            row1[8:16], row1[16:24], row1[24:32], row1[32:40], row1[40:48], row1[48:56], row1[56:64], row1[64:72]]

        card = BDFCard(card_fields, has_none=True)

        name = string(card, 1, 'name')
        #label_fields = [labeli for labeli in card[2:8] if labeli is not None]
        #label = ''.join(label_fields).strip()
        # assert len(label) <= 56, label

        axes = parse_components(card, 9, 'axes')
        comp = str(integer(card, 10, 'comp'))
        cp = integer_or_blank(card, 11, 'cp', default=0)
        xyz = [
            double_or_blank(card, 12, 'x', default=0.0),
            double_or_blank(card, 13, 'y', default=0.0),
            double_or_blank(card, 14, 'z', default=0.0),
        ]
        cd = integer_or_blank(card, 15, 'cd', default=cp)
        ind_dof = components_or_blank(card, 16, 'ind_dof', default='123')
        #return MONDSP1(name, label, axes, comp, xyz, cp=cp,
                       #cd=cd, ind_dof=ind_dof, comment=comment)
        self.cards.append((name, label, axes, comp, xyz, cp,
                           cd, ind_dof, comment))
        self.n += 1
        return self.n - 1

    @VectorizedBaseCard.parse_cards_check
    def parse_cards(self) -> None:
        ncards = len(self.cards)
        name = np.zeros(ncards, dtype='|U8')
        label = np.zeros(ncards, dtype='|U72')
        axes = np.zeros(ncards, dtype='int32')
        component = np.zeros(ncards, dtype='|U6')
        ind_dof = np.zeros(ncards, dtype='int32')
        cp = np.zeros(ncards, dtype='int32')
        xyz = np.zeros((ncards, 3), dtype='float64')
        cd = np.zeros(ncards, dtype='int32')
        for icard, card in enumerate(self.cards):
            (namei, labeli, axesi, compi, xyzi, cpi,
             cdi, ind_dofi, comment) = card
            name[icard] = namei
            label[icard] = labeli
            axes[icard] = axesi
            component[icard] = compi
            ind_dof[icard] = ind_dofi
            cp[icard] = cpi
            xyz[icard, :] = xyzi
            cd[icard] = cdi
        self._save(name, label, axes, component, ind_dof, cp, xyz, cd)
        #self.sort()
        self.cards = []

    def _save(self, name, label, axes, component, ind_dof, cp, xyz, cd):
        assert len(self.name) == 0, self.name
        self.name = name
        self.label = label
        self.axes = axes
        self.component = component
        self.ind_dof = ind_dof
        self.cp = cp
        self.xyz = xyz
        self.cd = cd

    def __apply_slice__(self, monitor: MONDSP1, i: np.ndarray) -> None:  # ignore[override]
        monitor.name = self.name[i]
        monitor.label = self.label[i]
        monitor.axes = self.axes[i]
        monitor.component = self.component[i]
        monitor.ind_dof = self.ind_dof[i]
        monitor.cp = self.cp[i]
        monitor.xyz = self.xyz[i, :]
        monitor.cd = self.cd[i]
        monitor.n = len(i)

    def set_used(self, used_dict: dict[str, np.ndarray]) -> None:
        used_dict['coord_id'].append(self.cp)
        used_dict['coord_id'].append(self.cd)

    def geom_check(self, missing: dict[str, np.ndarray]):
        #mids = hstack_msg([prop.material_id for prop in self.allowed_materials],
                          #msg=f'no materials for {self.type}')
        #mids.sort()
        coords = self.model.coord.coord_id
        #all_aecomp_names = self.model.aecomp.name
        #aecomp_names = np.unique(self.comp)
        ucoords = np.unique(np.hstack([self.cp, self.cd]))
        geom_check(self,
                   missing,
                   coord=(coords, ucoords),
                   #aecomp=(all_aecomp_names, aecomp_names),
                   )

    @property
    def max_id(self) -> int:
        return max(self.cp.max(), self.cd.max())

    def write_file(self, bdf_file: TextIOLike,
              size: int=8, is_double: bool=False,
              write_card_header: bool=False) -> None:
        if len(self.name) == 0:
            return None

        size = 8
        print_card, size = get_print_card_size(size, self.max_id)
        assert size == 8, size


        cps = array_default_int(self.cp, size=size)
        cds = array_str(self.cd, size=size)
        components = array_str(self.component, size=size)
        ind_dofs = array_str(self.ind_dof, size=size)
        is_default_cd = (self.cp == self.cd)
        cds[is_default_cd] = ''
        #cds = array_default_int(self.cd, default=0, size=size)
        xyzs = array_float(self.xyz, size=size, is_double=False).tolist()
        for name, label, axes, comp, ind_dof, \
            cp, xyz, cd in zip(self.name, self.label, self.axes,
                               components, ind_dofs,
                               cps, xyzs, cds):
            x, y, z = xyz

            msg = 'MONDSP1 %-8s%s\n' % (name, label)
            msg += '        %-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n' % (
                axes, comp, cp, x, y, z, cd, ind_dof)
            bdf_file.write(msg.rstrip('\n ') + '\n')
        return
