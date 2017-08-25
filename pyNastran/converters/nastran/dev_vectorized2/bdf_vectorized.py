from __future__ import print_function
from collections import defaultdict
from six import iteritems
import numpy as np
from pyNastran.bdf.bdf import BDF as BDF_, GRID

from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_blank, double_or_blank,
    components_or_blank)
from pyNastran.bdf.field_writer_8 import print_float_8, set_string8_blank_if_default
from pyNastran.bdf.field_writer_16 import print_float_16, set_string16_blank_if_default
from pyNastran.bdf.field_writer_double import print_scientific_double

from pyNastran.converters.nastran.dev_vectorized2.springs import (
    CELAS1, CELAS2, CELAS3, CELAS4, Springs)
from pyNastran.converters.nastran.dev_vectorized2.shells import CQUAD4v, CTRIA3v, Shells
from pyNastran.converters.nastran.dev_vectorized2.solids import (
    CTETRA4v, CPENTA6v, CHEXA8v, CPYRAM5v,
    CTETRA10v, CPENTA15v, CHEXA20v, CPYRAM13v, Solids)


class GRIDv(object):
    card_name = 'GRID'
    def __init__(self, model):
        self.model = model
        self._is_current = False
        self.nid = np.array([], dtype='int32')
        self.xyz = np.array([], dtype='float64')
        self.cp = np.array([], dtype='int32')
        self.cd = np.array([], dtype='int32')
        self.ps = np.array([], dtype='|U8')
        self.seid = np.array([], dtype='int32')

        self._nid = []
        self._xyz = []
        self._cp = []
        self._cd = []
        self._ps = []
        self._seid = []
        self.comment = defaultdict(str)

    def add(self, nid, xyz, cp=0, cd=0, ps='', seid=0, comment=''):
        # type: (int, Union[None, List[float], np.ndarray], int, int, str, int, str) -> None
        """
        Creates the GRID card

        Parameters
        ----------
        nid : int
            node id
        cp : int; default=0
            the xyz coordinate frame
        xyz : (3, ) float ndarray; default=None -> [0., 0., 0.]
            the xyz/r-theta-z/rho-theta-phi values
        cd : int; default=0
            the analysis coordinate frame
        ps : str; default=''
            Additional SPCs in the analysis coordinate frame (e.g. '123').
            This corresponds to DOF set ``SG``.
        seid : int; default=0
            superelement id
            TODO: how is this used by Nastran???
        comment : str; default=''
            a comment for the card
        """
        self.is_current = False
        self._nid.append(nid)
        self._xyz.append(xyz)
        self._cp.append(cp)
        self._cd.append(cd)
        self._ps.append(ps)
        self._seid.append(seid)
        if comment:
            self.comment[nid] = comment

    def add_card(self, card, comment=''):
        # type: (Any, str) -> GRID
        """
        Adds a GRID card from ``BDF.add_card(...)``

        Parameters
        ----------
        card : BDFCard()
            a BDFCard object
        comment : str; default=''
            a comment for the card
        """
        nfields = len(card)
        #print('card = %s' % card)
        #: Node ID
        nid = integer(card, 1, 'nid')

        #: Grid point coordinate system
        cp = integer_or_blank(card, 2, 'cp', 0)

        #: node location in local frame
        xyz = [
            double_or_blank(card, 3, 'x1', 0.),
            double_or_blank(card, 4, 'x2', 0.),
            double_or_blank(card, 5, 'x3', 0.)]

        if nfields > 6:
            #: Analysis coordinate system
            cd = integer_or_blank(card, 6, 'cd', 0)

            #: SPC constraint
            ps = components_or_blank(card, 7, 'ps', '')
            #u(integer_or_blank(card, 7, 'ps', ''))

            #: Superelement ID
            seid = integer_or_blank(card, 8, 'seid', 0)
            assert len(card) <= 9, 'len(GRID card) = %i\ncard=%s' % (len(card), card)
        else:
            cd = 0
            ps = ''
            seid = 0
        self.add(nid, xyz, cp, cd, ps, seid, comment=comment)

    def check_if_current(self, nid, nids):
        """we split this up to reason about it easier"""
        if self._is_current:
            if nid in nids:
                # card exists, so we use that slot
                add_card = False
            else:
                add_card = True
        else:
            add_card = True
        return add_card

    def update(self, grid):
        """functions like a dictionary"""
        nid = grid.nid

        add_grid = self.check_if_current(nid, self.nid)
        if add_grid:
            self.add(nid, grid.xyz, cp=grid.cp, cd=grid.cd,
                     ps=grid.ps, seid=grid.seid, comment=grid.comment)
            self._is_current = False
        else:
            inid = np.where(nid == self.nid)[0]
            self.nid[inid] = grid.nid
            self.xyz[inid] = grid.xyz
            self.cp[inid] = grid.cp
            self.cd[inid] = grid.cd
            self.ps[inid] = grid.ps
            self.seid[inid] = grid.seid
            #self.comment[nid] = comment
            #self._is_current = True  # implicit

    def _make_current(self):
        """creates an array of the GRID points"""
        if not self._is_current:
            if len(self.nid) > 0: # there are already nodes in self.nid
                self.nid = np.hstack([self.nid, self._nid])
                self.xyz = np.vstack([self.xyz, self._xyz])
                self.cp = np.hstack([self.cp, self._cp])
                self.cd = np.hstack([self.cd, self._cd])
                self.ps = np.hstack([self.ps, self._ps])
                self.seid = np.hstack([self.seid, self._seid])
                # don't need to handle comments
            else:
                self.nid = np.array(self._nid)
                self.xyz = np.array(self._xyz)
                self.cp = np.array(self._cp)
                self.cd = np.array(self._cd)
                self.ps = np.array(self._ps)
                self.seid = np.array(self._seid)
            assert len(self.nid) == len(np.unique(self.nid))
            #print(self.nid)
            self._nid = []
            self._xyz = []
            self._cp = []
            self._cd = []
            self._ps = []
            self._seid = []
            self._is_current = True
        #else:
            #print('no GRIDs')

    def cross_reference(self, model):
        """does this do anything?"""
        self._make_current()

    def write_card(self, size=8, is_double=False):
        # type: (int, bool) -> str
        """
        The writer method used by BDF.write_card

        Parameters
        ----------
        size : int; default=8
            the size of the card (8/16)
        is_double : bool; default=False
            should this card be written with double precision

        Returns
        -------
        msg : str
            the card as a string
        """
        if size == 8:
            return self.write_card_8()
        else:
            return self.write_card_16(is_double)

    def write_card_8(self):
        # type: () -> str
        """
        Writes a GRID card in 8-field format
        """
        self._make_current()
        msg = ''
        for nid, xyz, cp, cd, ps, seid in zip(self.nid, self.xyz, self.cp, self.cd, self.ps, self.seid):
            cps = set_string8_blank_if_default(cp, 0)
            if [cd, ps, seid] == [0, '', 0]:
                # default
                print_float_8(xyz[0])
                print_float_8(xyz[1])
                print_float_8(xyz[2])
                msgi = 'GRID    %8i%8s%s%s%s\n' % (
                    nid, cps,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                )
                msg += self.comment[nid] + msgi
            else:
                cds = set_string8_blank_if_default(cd, 0)
                seids = set_string8_blank_if_default(seid, 0)
                msgi = 'GRID    %8i%8s%s%s%s%s%8s%s\n' % (
                    nid, cps,
                    print_float_8(xyz[0]),
                    print_float_8(xyz[1]),
                    print_float_8(xyz[2]),
                    cds, ps, seids)
                msg += self.comment[nid] + msgi
        return msg

    def write_card_16(self, is_double=False):
        # type: (bool) -> str
        """
        Writes a GRID card in 16-field format
        """
        self._make_current()
        xyz = self.xyz
        cp = set_string16_blank_if_default(cp, 0)
        cd = set_string16_blank_if_default(cd, 0)
        seid = set_string16_blank_if_default(seid, 0)

        if is_double:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_scientific_double(xyz[0]),
                           print_scientific_double(xyz[1]),
                           print_scientific_double(xyz[2]),
                           cd, self.ps, seid))
        else:
            if [cd, self.ps, self.seid] == [0, '', 0]:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2])))
            else:
                msg = ('GRID*   %16i%16s%16s%16s\n'
                       '*       %16s%16s%16s%16s\n' % (
                           self.nid,
                           cp,
                           print_float_16(xyz[0]),
                           print_float_16(xyz[1]),
                           print_float_16(xyz[2]),
                           cd, self.ps, seid))
        return self.comment + msg

    def __repr__(self):
        self._make_current()
        msg = 'GRID_Vector; ngrids=%s:\n' % len(self.nid)
        msg += '  nid = %s\n' % self.nid

        ucp = np.unique(self.cp)
        if len(ucp) == 1 and ucp[0] == 0:
            msg += '  ucp = %s\n' % ucp
        else:
            msg += '  cp = %s\n' % self.cp

        ucd = np.unique(self.cd)
        if len(ucd) == 1 and ucd[0] == 0:
            msg += '  ucd = %s\n' % ucd
        else:
            msg += '  cd = %s\n' % self.cd

        ups = np.unique(self.ps)
        if len(ups) == 1 and ups[0] == '':
            msg += '  ups = %s\n' % ups
        else:
            msg += '  ps = %s\n' % self.ps

        useid = np.unique(self.seid)
        if len(useid) == 1 and useid[0] == 0:
            msg += '  useid = %s\n' % useid
        else:
            msg += '  seid = %s\n' % self.seid
        #msg += '  xyz =\n%s' % self.xyz
        return msg
    #def __iter__(self):
        #pass
    #def __next__(self):
        #pass
    #def __items__(self):
        #pass
    #def __keys__(self):
        #pass
    #def __values__(self):
        #pass
    def __getitem__(self, i):
        """this works on index"""
        self._make_current()
        nid = self.nid[i]
        return GRID(nid, self.xyz[i], cp=self.cp[i], cd=self.cd[i],
                    ps=self.ps[i], seid=self.seid[i], comment=self.comment[nid])

    def get_grid_by_nid(self, nid):
        self._make_current()
        inid = np.searchsorted(self.nid, nid)
        return self[inid]

    def __setitem__(self, i, value):
        pass
    #def __delitem__(self, i):
        #pass

class Elements(object):
    """stores all the elements"""
    def __init__(self, model):
        self.model = model
        self.shells = model.shells
        self.solids = model.solids
        self.springs = model.springs

    def repr_indent(self, indent=''):
        indent = ''
        msg = '%s<Elements> : nelements=%s\n' % (indent, len(self))

        msg += '%s  springs:  %s\n' % (indent, len(self.springs))
        msg += self.springs.repr_indent('    ')

        msg += '%s  shells:  %s\n' % (indent, len(self.shells))
        msg += self.shells.repr_indent('    ')

        msg += '%s  solids:  %s\n' % (indent, len(self.solids))
        msg += self.solids.repr_indent('    ')
        return msg

    def write_card(self, size=8, is_double=False, bdf_file=None):
        assert bdf_file is not None
        self.springs.write_card(size, is_double, bdf_file)
        self.shells.write_card(size, is_double, bdf_file)
        self.solids.write_card(size, is_double, bdf_file)

    def __len__(self):
        return len(self.springs) + len(self.shells) + len(self.solids)
    def __repr__(self):
        return self.repr_indent('')


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

        self.ctria3 = CTRIA3v(model)
        self.cquad4 = CQUAD4v(model)
        self.ctria6 = CTRIA3v(model)  # TODO: temp
        self.cquad8 = CQUAD4v(model)  # TODO: temp
        #self.cquad = CQUADv(model)  # TODO: temp
        self.shells = Shells(model)
        #self.pshell = PSHELLv(model)  # TODO: temp

        #self.ctriax = CTRIA3v(model)   # TODO: temp
        #self.cquadx = CTRIA3v(model)   # TODO: temp
        #self.ctriax6 = CTRIA3v(model)  # TODO: temp
        #self.cquadx8 = CTRIA3v(model)  # TODO: temp

        #self.crod = CRODv(model)      # TODO: temp
        #self.conrod = CONRODv(model)  # TODO: temp
        #self.ctube = CTUBEv(model)    # TODO: temp
        #self.rods = Rods(model)       # TODO: temp

        self.celas1 = CELAS1(model)
        self.celas2 = CELAS2(model)
        self.celas3 = CELAS3(model)
        self.celas4 = CELAS4(model)
        self.springs = Springs(model)

        #self.cdamp1 = CDAMP1(model)    # TODO: temp
        #self.cdamp2 = CDAMP2(model)    # TODO: temp
        #self.cdamp3 = CDAMP3(model)    # TODO: temp
        #self.cdamp4 = CDAMP4(model)    # TODO: temp
        #self.cdamp5 = CDAMP5(model)    # TODO: temp
        #self.dampers = Dampers(model)  # TODO: temp

        self.ctetra4 = CTETRA4v(model)
        self.ctetra10 = CTETRA10v(model)
        self.chexa8 = CHEXA8v(model)
        self.chexa20 = CHEXA20v(model)
        self.cpenta6 = CPENTA6v(model)
        self.cpenta15 = CPENTA15v(model)
        self.cpyram5 = CPYRAM5v(model)
        self.cpyram13 = CPYRAM13v(model)
        self.solids = Solids(model)

        self.elements2 = Elements(model)

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

    def _prepare_cquad4(self, card, card_obj, comment=''):
        self.cquad4.add_card(card_obj, comment=comment)
    def _prepare_ctria3(self, card, card_obj, comment=''):
        self.ctria3.add_card(card_obj, comment=comment)

    def _prepare_celas1(self, card, card_obj, comment=''):
        self.celas1.add_card(card_obj, comment=comment)
    def _prepare_celas2(self, card, card_obj, comment=''):
        self.celas2.add_card(card_obj, comment=comment)
    def _prepare_celas3(self, card, card_obj, comment=''):
        self.celas3.add_card(card_obj, comment=comment)
    def _prepare_celas4(self, card, card_obj, comment=''):
        self.celas4.add_card(card_obj, comment=comment)

    #def _prepare_ctria6(self, card, card_obj, comment=''):
        #self.ctria6.add_card(card_obj, comment=comment)
    #def _prepare_cquad8(self, card, card_obj, comment=''):
        #self.cquad8.add_card(card_obj, comment=comment)
    #def _prepare_ctetra(self, card, card_obj, comment=''):
        #self.ctetra.add_card(card_obj, comment=comment)

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

    def _update_card_parser(self):
        del self._card_parser['GRID']
        self._card_parser_prepare['GRID'] = self._prepare_grid

        del self._card_parser['CQUAD4']
        self._card_parser_prepare['CQUAD4'] = self._prepare_cquad4

        del self._card_parser['CTRIA3']
        self._card_parser_prepare['CTRIA3'] = self._prepare_ctria3

        del self._card_parser['CELAS1']
        del self._card_parser['CELAS2']
        del self._card_parser['CELAS3']
        del self._card_parser['CELAS4']
        self._card_parser_prepare['CELAS1'] = self._prepare_celas1
        self._card_parser_prepare['CELAS2'] = self._prepare_celas2
        self._card_parser_prepare['CELAS3'] = self._prepare_celas3
        self._card_parser_prepare['CELAS4'] = self._prepare_celas4

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
        raise NotImplementedError('grids')
        #return self.grid

    @grids.setter
    def grids(self, grid):
        print('setter')
        key = grid.nid
        self.grid.update(grid)
        #self.grid.add_grid(grid.nid, grid.xyz, cp=grid.cp, cd=grid.cd,
                           #ps=grid.ps, seid=grid.seid, comment=grid.comment)

    def _write_nodes(self, bdf_file, size=8, is_double=False):
        # type: (Any, int, bool) -> None
        """
        Writes the NODE-type cards
        """
        BDF_._write_nodes(self, bdf_file, size=size, is_double=is_double)
        bdf_file.write(self.grid.write_card(size=size, is_double=is_double))

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
    return model


def test_vectorized():
    import os
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    model_path = os.path.join(pkg_path, '../', 'models')

    bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
    bdf_filename = os.path.join(model_path, 'bwb', 'BWB_saero.bdf')
    bdf_filename = os.path.join(model_path, 'elements', 'static_elements.bdf')
    model = read_bdf(bdf_filename, validate=True, xref=False, punch=False,
                     skip_cards=None, encoding=None, log=None, debug=True,
                     mode='msc')

    #model.grids[10] = GRID(10, [0., 0., 0.])
    print(model.grid)
    print(model.grid[10]) # nid or index?
    print(model.grid.get_grid_by_nid(10))
    print(model.cquad4)
    print(model.ctria3)
    print(model.shells)
    print(model.solids)
    print(model.ctetra4)
    print(model.ctetra10)

    print(model.celas1)
    print(model.celas2)
    print(model.celas3)
    print(model.celas4)

    print(model.elements2)
    print(len(model.elements2))
    out_filename = 'spike.bdf'
    model.write_bdf(out_filename, encoding=None, size=8, is_double=False,
                    interspersed=False, enddata=None,
                    close=True)

if __name__ == '__main__':  # pragma: no cover
    test_vectorized()
