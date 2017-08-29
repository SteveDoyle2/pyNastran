from __future__ import print_function
from collections import defaultdict
from six import iteritems
import numpy as np

from pyNastran.bdf.bdf import BDF as BDF_
from pyNastran.converters.nastran.dev_vectorized2.cards.nodes import GRIDv, Nodes
from pyNastran.converters.nastran.dev_vectorized2.cards.elements import Elements

from pyNastran.converters.nastran.dev_vectorized2.cards.springs import (
    CELAS1, CELAS2, CELAS3, CELAS4, Springs)
from pyNastran.converters.nastran.dev_vectorized2.cards.dampers import (
    CDAMP1, CDAMP2, CDAMP3, CDAMP4, Dampers)
from pyNastran.converters.nastran.dev_vectorized2.cards.bush import (
    CBUSHv, Bushes)

from pyNastran.converters.nastran.dev_vectorized2.cards.rods import (
    CONRODv, CRODv, CTUBEv, Rods)
from pyNastran.converters.nastran.dev_vectorized2.cards.loads import (
    Loads, PLOAD2v, PLOAD4v, FORCEv, FORCE1v, FORCE2v)
from pyNastran.converters.nastran.dev_vectorized2.cards.bars import CBARv, Bars
from pyNastran.converters.nastran.dev_vectorized2.cards.beams import CBEAMv, Beams
from pyNastran.converters.nastran.dev_vectorized2.cards.shears import CSHEARv, Shears
from pyNastran.converters.nastran.dev_vectorized2.cards.shells import (
    CTRIA3v, CTRIA6v, CQUAD4v, CQUAD8v, Shells)
from pyNastran.converters.nastran.dev_vectorized2.cards.solids import (
    CTETRA4v, CPENTA6v, CHEXA8v, CPYRAM5v,
    CTETRA10v, CPENTA15v, CHEXA20v, CPYRAM13v, Solids)


class BDF(BDF_):
    """
    Uses the BDF class, but overwrites a few classes/methods

    Vectorized:
     - GRID
    """
    def __init__(self, debug=True, log=None, mode='msc'):
        # type: (Optional[bool], SimpleLogger, str) -> None
        BDF_.__init__(self, debug=debug, log=log, mode=mode)
        #super(BDF, self).__init__(debug=debug, log=log, mode=mode)
        #del
        #self._grids_temp = []

        model = self
        self.grid = GRIDv(model)
        self.nodes2 = Nodes(model)

        self.celas1 = CELAS1(model)
        self.celas2 = CELAS2(model)
        self.celas3 = CELAS3(model)
        self.celas4 = CELAS4(model)
        self.springs = Springs(model)

        self.cdamp1 = CDAMP1(model)
        self.cdamp2 = CDAMP2(model)
        self.cdamp3 = CDAMP3(model)
        self.cdamp4 = CDAMP4(model)
        #self.cdamp5 = CDAMP5(model)    # TODO: temp
        self.dampers = Dampers(model)

        self.cbush = CBUSHv(model)
        self.bushes = Bushes(model)

        self.crod = CRODv(model)
        self.conrod = CONRODv(model)
        self.ctube = CTUBEv(model)
        self.rods = Rods(model)

        self.cbar = CBARv(model)
        self.bars = Bars(model)

        self.cbeam = CBEAMv(model)  # TODO: temp
        self.beams = Beams(model)   # TODO: temp

        self.ctria3 = CTRIA3v(model)
        self.cquad4 = CQUAD4v(model)
        self.ctria6 = CTRIA6v(model)
        self.cquad8 = CQUAD8v(model)
        self.cquad = CQUAD4v(model)   # TODO: temp
        self.shells = Shells(model)
        #self.pshell = PSHELLv(model)  # TODO: temp

        self.cshear = CSHEARv(model)
        self.shears = Shears(model)

        #self.ctriax = CTRIA3v(model)   # TODO: temp
        #self.cquadx = CTRIA3v(model)   # TODO: temp
        #self.ctriax6 = CTRIA3v(model)  # TODO: temp
        #self.cquadx8 = CTRIA3v(model)  # TODO: temp

        self.ctetra4 = CTETRA4v(model)
        self.ctetra10 = CTETRA10v(model)
        self.chexa8 = CHEXA8v(model)
        self.chexa20 = CHEXA20v(model)
        self.cpenta6 = CPENTA6v(model)
        self.cpenta15 = CPENTA15v(model)
        self.cpyram5 = CPYRAM5v(model)
        self.cpyram13 = CPYRAM13v(model)
        self.solids = Solids(model)

        self.elements2 = Elements(model)  # TODO: change this name


        self.force = FORCEv(model)
        self.force1 = FORCE1v(model)
        self.force2 = FORCE2v(model)
        self.pload2 = PLOAD2v(model)
        self.pload4 = PLOAD4v(model)
        self.loads2 = Loads(model)

        self._update_card_parser()

    def _prepare_grid(self, card, card_obj, comment=''):
        self.grid.add_card(card_obj, comment=comment)
    def _add_node_object(self, node, allow_overwrites=False):
        raise AttributeError("'BDF' object has no attribute '_add_node_object'")
    def _add_element_object(self, elem, allow_overwrites=False):
        key = elem.eid
        assert key > 0, 'eid=%s must be positive; elem=\n%s' % (key, elem)
        #if key in self.shells.eids:
        if key in self.elements and not allow_overwrites:
            if not elem == self.elements[key]:
                self._duplicate_elements.append(elem)
                if self._stop_on_duplicate_error:
                    self.pop_parse_errors()
        else:
            self.elements[key] = elem
            self._type_to_id_map[elem.type].append(key)

    def _prepare_celas1(self, card, card_obj, comment=''):
        self.celas1.add_card(card_obj, comment=comment)
    def _prepare_celas2(self, card, card_obj, comment=''):
        self.celas2.add_card(card_obj, comment=comment)
    def _prepare_celas3(self, card, card_obj, comment=''):
        self.celas3.add_card(card_obj, comment=comment)
    def _prepare_celas4(self, card, card_obj, comment=''):
        self.celas4.add_card(card_obj, comment=comment)

    def _prepare_cdamp1(self, card, card_obj, comment=''):
        self.cdamp1.add_card(card_obj, comment=comment)
    def _prepare_cdamp2(self, card, card_obj, comment=''):
        self.cdamp2.add_card(card_obj, comment=comment)
    def _prepare_cdamp3(self, card, card_obj, comment=''):
        self.cdamp3.add_card(card_obj, comment=comment)
    def _prepare_cdamp4(self, card, card_obj, comment=''):
        self.cdamp4.add_card(card_obj, comment=comment)
    #def _prepare_cdamp5(self, card, card_obj, comment=''):
        #self.cdamp5.add_card(card_obj, comment=comment)

    def _prepare_cbush(self, card, card_obj, comment=''):
        self.cbush.add_card(card_obj, comment=comment)
    #def _prepare_cbush1d(self, card, card_obj, comment=''):
        #self.cbush1d.add_card(card_obj, comment=comment)
    #def _prepare_cbush2d(self, card, card_obj, comment=''):
        #self.cbush2d.add_card(card_obj, comment=comment)

    def _prepare_conrod(self, card, card_obj, comment=''):
        self.conrod.add_card(card_obj, comment=comment)
    def _prepare_crod(self, card, card_obj, comment=''):
        self.crod.add_card(card_obj, comment=comment)
    def _prepare_ctube(self, card, card_obj, comment=''):
        self.ctube.add_card(card_obj, comment=comment)

    def _prepare_cbar(self, card, card_obj, comment=''):
        self.cbar.add_card(card_obj, comment=comment)
    def _prepare_cbeam(self, card, card_obj, comment=''):
        self.cbeam.add_card(card_obj, comment=comment)

    def _prepare_cquad4(self, card, card_obj, comment=''):
        self.cquad4.add_card(card_obj, comment=comment)
    def _prepare_ctria3(self, card, card_obj, comment=''):
        self.ctria3.add_card(card_obj, comment=comment)
    def _prepare_ctria6(self, card, card_obj, comment=''):
        self.ctria6.add_card(card_obj, comment=comment)
    def _prepare_cquad8(self, card, card_obj, comment=''):
        self.cquad8.add_card(card_obj, comment=comment)
    #def _prepare_cquad(self, card, card_obj, comment=''):
        #self.cquad.add_card(card_obj, comment=comment)

    def _prepare_cshear(self, card, card_obj, comment=''):
        self.cshear.add_card(card_obj, comment=comment)

    def _prepare_ctetra(self, card, card_obj, comment=''):
        """adds a CTETRA4/CTETRA10"""
        if len(card_obj) == 7:
            self.ctetra4.add_card(card_obj, comment=comment)
        else:
            self.ctetra10.add_card(card_obj, comment=comment)

    def _prepare_cpenta(self, card, card_obj, comment=''):
        """adds a CPENTA6/CPENTA15"""
        if len(card_obj) == 9:
            self.cpenta6.add_card(card_obj, comment=comment)
        else:
            self.cpenta15.add_card(card_obj, comment=comment)

    def _prepare_chexa(self, card, card_obj, comment=''):
        """adds a CHEXA8/CHEXA20"""
        if len(card_obj) == 11:
            self.chexa8.add_card(card_obj, comment=comment)
        else:
            self.chexa20.add_card(card_obj, comment=comment)

    def _prepare_cpyram(self, card, card_obj, comment=''):
        """adds a CPYRAM5/CPYRAM13"""
        if len(card_obj) == 8:
            self.cpyram5.add_card(card_obj, comment=comment)
        else:
            self.cpyram13.add_card(card_obj, comment=comment)

    def _prepare_force(self, card, card_obj, comment=''):
        self.force.add_card(card_obj, comment=comment)
    def _prepare_force1(self, card, card_obj, comment=''):
        self.force1.add_card(card_obj, comment=comment)
    def _prepare_force2(self, card, card_obj, comment=''):
        self.force2.add_card(card_obj, comment=comment)
    def _prepare_pload2(self, card, card_obj, comment=''):
        self.pload2.add_card(card_obj, comment=comment)
    def _prepare_pload4(self, card, card_obj, comment=''):
        self.pload4.add_card(card_obj, comment=comment)

    def _update_card_parser(self):
        del self._card_parser['GRID']
        self._card_parser_prepare['GRID'] = self._prepare_grid

        del self._card_parser['CELAS1']
        del self._card_parser['CELAS2']
        del self._card_parser['CELAS3']
        del self._card_parser['CELAS4']
        self._card_parser_prepare['CELAS1'] = self._prepare_celas1
        self._card_parser_prepare['CELAS2'] = self._prepare_celas2
        self._card_parser_prepare['CELAS3'] = self._prepare_celas3
        self._card_parser_prepare['CELAS4'] = self._prepare_celas4

        del self._card_parser['CDAMP1']
        del self._card_parser['CDAMP2']
        del self._card_parser['CDAMP3']
        #del self._card_parser['CDAMP4']  # no
        self._card_parser_prepare['CDAMP1'] = self._prepare_cdamp1
        self._card_parser_prepare['CDAMP2'] = self._prepare_cdamp2
        self._card_parser_prepare['CDAMP3'] = self._prepare_cdamp3
        self._card_parser_prepare['CDAMP4'] = self._prepare_cdamp4

        del self._card_parser['CBUSH']
        #del self._card_parser['CBUSH1D']
        #del self._card_parser['CBUSH2D']
        self._card_parser_prepare['CBUSH'] = self._prepare_cbush
        #self._card_parser_prepare['CBUSH1D'] = self._prepare_cbush1d
        #self._card_parser_prepare['CBUSH2D'] = self._prepare_cbush2d

        del self._card_parser['CONROD']
        del self._card_parser['CROD']
        del self._card_parser['CTUBE']
        self._card_parser_prepare['CONROD'] = self._prepare_conrod
        self._card_parser_prepare['CROD'] = self._prepare_crod
        self._card_parser_prepare['CTUBE'] = self._prepare_ctube

        del self._card_parser['CBAR']
        del self._card_parser['CBEAM']
        self._card_parser_prepare['CBAR'] = self._prepare_cbar
        self._card_parser_prepare['CBEAM'] = self._prepare_cbeam

        del self._card_parser['CTRIA3']
        del self._card_parser['CTRIA6']
        del self._card_parser['CQUAD4']
        del self._card_parser['CQUAD8']
        #del self._card_parser['CQUAD']
        self._card_parser_prepare['CTRIA3'] = self._prepare_ctria3
        self._card_parser_prepare['CTRIA6'] = self._prepare_ctria6
        self._card_parser_prepare['CQUAD4'] = self._prepare_cquad4
        self._card_parser_prepare['CQUAD8'] = self._prepare_cquad8
        #self._card_parser_prepare['CQUAD'] = self._prepare_cquad

        del self._card_parser['CSHEAR']
        self._card_parser_prepare['CSHEAR'] = self._prepare_cshear


        del self._card_parser['PLOAD2']
        del self._card_parser['PLOAD4']
        del self._card_parser['FORCE']
        del self._card_parser['FORCE1']
        del self._card_parser['FORCE2']
        self._card_parser_prepare['PLOAD2'] = self._prepare_pload2
        self._card_parser_prepare['PLOAD4'] = self._prepare_pload4
        self._card_parser_prepare['FORCE'] = self._prepare_force
        self._card_parser_prepare['FORCE1'] = self._prepare_force1
        self._card_parser_prepare['FORCE2'] = self._prepare_force2


    #def add_grid(self, nid, xyz, cp=0, cd=0, ps='', seid=0, comment=''):
        #pass
    #def cross_reference(self, xref=True, xref_nodes=True, xref_elements=True,
                       #xref_nodes_with_elements=True,
                       #xref_properties=True,
                       #xref_masses=True,
                       #xref_materials=True,
                       #xref_loads=True,
                       #xref_constraints=True,
                       #xref_aero=True, xref_sets=True,
                       #xref_optimization=True):
        #self.grid.cross_reference(self)

#def _add_node_object(self, node, allow_overwrites=False):
    ## type: (Any, bool) -> None
    #"""adds a GRID card"""
    #key = node.nid
    #if key in self.nodes and not allow_overwrites:
        #if not node == self.nodes[key]:
            #assert node.nid not in self.nodes, 'nid=%s\nold_node=\n%snew_node=\n%s' % (node.nid, self.nodes[key], node)
        #else:
            ##print('GRID was duplicated...nid=%s; node=\n%s' % (key, node))
            #pass
    #else:
        #assert key > 0, 'nid=%s node=%s' % (key, node)
        #self.node.add_[key] = node
        #self._type_to_id_map[node.type].append(key)

    @property
    def grids(self):
        raise NotImplementedError('grids is removed; use grid')
        #return self.grid

    @grids.setter
    def grids(self, grid):
        print('setter')
        key = grid.nid
        self.grid.update(grid)
        #self.grid.add_grid(grid.nid, grid.xyz, cp=grid.cp, cd=grid.cd,
                           #ps=grid.ps, seid=grid.seid, comment=grid.comment)

    def xyz_cid0(self):
        return self.xyz_cid0

    def _write_nodes(self, bdf_file, size=8, is_double=False):
        # type: (Any, int, bool) -> None
        """
        Writes the NODE-type cards
        """
        BDF_._write_nodes(self, bdf_file, size=size, is_double=is_double)
        self.nodes2.write_card(size=size, is_double=is_double, bdf_file=bdf_file)

    def _write_elements(self, bdf_file, size=8, is_double=False):
        # type: (Any, int, bool) -> None
        """
        Writes the elements in a sorted order
        """
        if self.elements:
            bdf_file.write('$ELEMENTS\n')
            if self.is_long_ids:
                for (eid, element) in sorted(iteritems(self.elements)):
                    bdf_file.write(element.write_card_16(is_double))
            else:
                for (eid, element) in sorted(iteritems(self.elements)):
                    try:
                        bdf_file.write(element.write_card(size, is_double))
                    except:
                        print('failed printing element...'
                              'type=%s eid=%s' % (element.type, eid))
                        raise
        self.elements2.write_card(size, is_double, bdf_file)
        #bdf_file.write(self.shells.write_card(size, is_double))
        #bdf_file.write(self.solids.write_card(size, is_double))
        if self.ao_element_flags:
            for (eid, element) in sorted(iteritems(self.ao_element_flags)):
                bdf_file.write(element.write_card(size, is_double))
        self._write_nsm(bdf_file, size, is_double)

    def get_displacement_index_xyz_cp_cd(self, fdtype='float64', idtype='int32'):
        # type: (str, str, bool) -> Any
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        fdtype : str
            the type of xyz_cp
        idtype : str
            the type of nid_cp_cd

        Returns
        -------
        icd_transform : dict{int cd : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        nid_cp_cd : (n, 3) int ndarray
            node id, CP, CD for each node

        Example
        --------
        # assume GRID 1 has a CD=10, CP=0
        # assume GRID 2 has a CD=10, CP=0
        # assume GRID 5 has a CD=50, CP=0
        >>> model.point_ids
        [1, 2, 5]
        >>> out = model.get_displacement_index_xyz_cp_cd()
        >>> icd_transform, icp_transform, xyz_cp, nid_cp_cd = out
        >>> nid_cp_cd
        [
           [1, 0, 10],
           [2, 0, 10],
           [5, 0, 50],
        ]
        >>> icd_transform[10]
        [0, 1]

        >>> icd_transform[50]
        [2]
        """
        return self.nodes2.get_displacement_index_xyz_cp_cd(
            fdtype=fdtype, idtype=idtype)

    def get_node_index(self, nids, allow0=False):
        """maps the requested nodes to their desired index in the array"""
        i = self.nodes2.get_node_index(nids, allow0=allow0)
        return i

    def validate(self):
        pass
    def get_bdf_stats(self, return_type='string'):
        pass

def read_bdf(bdf_filename=None, validate=True, xref=True, punch=False,
             skip_cards=None,
             encoding=None, log=None, debug=True, mode='msc'):
    # type: (Union[str, None], bool, bool, bool, Union[List[str], None], Union[str, None], Union[SimpleLogger, None], Optional[bool], str) -> BDF
    """
    Creates the BDF object

    Parameters
    ----------
    bdf_filename : str (default=None -> popup)
        the bdf filename
    debug : bool/None
        used to set the logger if no logger is passed in
            True:  logs debug/info/error messages
            False: logs info/error messages
            None:  logs error messages
    log : logging module object / None
        if log is set, debug is ignored and uses the
        settings the logging object has
    validate : bool; default=True
        runs various checks on the BDF
    xref :  bool; default=True
        should the bdf be cross referenced
    punch : bool; default=False
        indicates whether the file is a punch file
    skip_cards : List[str]; default=None
        None : include all cards
        list of cards to skip
    encoding : str; default=None -> system default
        the unicode encoding
    mode : str; default='msc'
        the type of Nastran
        valid_modes = {'msc', 'nx'}

    Returns
    -------
    model : BDF()
        an BDF object

    .. code-block:: python

       >>> bdf = BDF()
       >>> bdf.read_bdf(bdf_filename, xref=True)
       >>> g1 = bdf.Node(1)
       >>> print(g1.get_position())
       [10.0, 12.0, 42.0]
       >>> bdf.write_card(bdf_filename2)
       >>> print(bdf.card_stats())

       ---BDF Statistics---
       SOL 101
       bdf.nodes = 20
       bdf.elements = 10
       etc.

    .. note :: this method will change in order to return an object that
               does not have so many methods
    .. todo:: finish this
    """
    model = BDF(log=log, debug=debug, mode=mode)
    if skip_cards:
        model.disable_cards(skip_cards)
    model.read_bdf(bdf_filename=bdf_filename, validate=validate,
                   xref=xref, punch=punch, read_includes=True, encoding=encoding)
    model.elements2.make_current()
    return model
