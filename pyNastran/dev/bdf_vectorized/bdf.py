# coding: utf-8
# pylint: disable=W0201,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF

see https://docs.plm.automation.siemens.com/tdoc/nxnastran/10/help/#uid:index
"""
import sys
import traceback
from pickle import load, dump
from collections import defaultdict
from typing import Optional

import numpy as np
from cpylog import SimpleLogger, get_logger, __version__ as CPYLOG_VERSION
assert CPYLOG_VERSION >= '1.6.0', CPYLOG_VERSION

from pyNastran.utils import object_attributes, check_path, PathLike
from pyNastran.bdf.bdf_interface.utils import (
    to_fields, _parse_pynastran_header, parse_executive_control_deck)

from pyNastran.bdf.utils import parse_patran_syntax

from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16

from pyNastran.bdf.cards.base_card import _format_comment
from pyNastran.bdf.cards.utils import wipe_empty_fields

#from pyNastran.bdf.write_path import write_include
from pyNastran.bdf.bdf_interface.assign_type import (
    integer, integer_or_string, string)

#from pyNastran.bdf.errors import CrossReferenceError, DuplicateIDsError, CardParseSyntaxError
#from pyNastran.bdf.field_writer_16 import print_field_16

from pyNastran.bdf.case_control_deck import CaseControlDeck

from pyNastran.bdf.bdf_interface.utils import fill_dmigs, _prep_comment
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.dev.bdf_vectorized.bdf_interface2.write_mesh import WriteMesh
from pyNastran.dev.bdf_vectorized.bdf_interface2.get_card import GetMethods
from pyNastran.dev.bdf_vectorized.bdf_interface2.cross_reference import CrossReference
from pyNastran.dev.bdf_vectorized.bdf_interface2.add_card import AddCard
from pyNastran.bdf.field_writer_16 import print_field_16

from pyNastran.dev.bdf_vectorized.cards.constraints.spc import SPC, get_spc_constraint
from pyNastran.dev.bdf_vectorized.cards.constraints.spcd import SPCD
from pyNastran.dev.bdf_vectorized.cards.constraints.spc1 import SPC1, get_spc1_constraint
from pyNastran.dev.bdf_vectorized.cards.constraints.spcadd import SPCADD, get_spcadd_constraint

from pyNastran.dev.bdf_vectorized.cards.constraints.mpc import MPC, get_mpc_constraint
#from pyNastran.dev.bdf_vectorized.cards.constraints.mpcax import MPCAX
from pyNastran.dev.bdf_vectorized.cards.constraints.mpcadd import MPCADD

from pyNastran.dev.bdf_vectorized.cards.deqatn import DEQATN
from pyNastran.dev.bdf_vectorized.cards.dynamic import (
    #DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ4,
    TSTEP, TSTEPNL, NLPARM, NLPCI, #TF
)

from pyNastran.dev.bdf_vectorized.cards.aero.aero_cards import (
    AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESTAT,
    AESURF, AESURFS, AERO, AEROS, CSSCHD,
    CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
    PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
    MONPNT1, MONPNT2, MONPNT3,
    FLFACT, FLUTTER, GUST, MKAERO1,
    MKAERO2, SPLINE1, SPLINE2, SPLINE3, SPLINE4,
    SPLINE5, TRIM, DIVERG)
from pyNastran.dev.bdf_vectorized.cards.optimization import (
    DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
    DRESP1, DRESP2, DRESP3,
    DVCREL1, DVCREL2,
    DVMREL1, DVMREL2,
    DVPREL1, DVPREL2,
    DVGRID)
from pyNastran.dev.bdf_vectorized.cards.bdf_sets import (
    ASET, BSET, CSET, QSET, USET,
    ASET1, BSET1, CSET1, QSET1, USET1,
    SET1, SET3, #RADSET,
    SEBSET, SECSET, SEQSET, # SEUSET
    SEBSET1, SECSET1, SEQSET1, # SEUSET1
    SESET, #SEQSEP,
)



# old cards
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE
from pyNastran.bdf.cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA
from pyNastran.bdf.cards.elements.elements import PLOTEL #CFAST, CGAP, CRAC2D, CRAC3D,
from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL
#from pyNastran.bdf.cards.loads.loads import (
    #DAREA, #LSEQ, SLOAD, DAREA, RANDPS, RFORCE, RFORCE1, SPCD, LOADCYN
#)
from pyNastran.bdf.errors import DuplicateIDsError, CrossReferenceError, CardParseSyntaxError
#from pyNastran.bdf.errors import (CrossReferenceError, DuplicateIDsError,
                                  #CardParseSyntaxError, UnsupportedCard, DisabledCardError,
                                  #SuperelementFlagError, ReplicationError)

from pyNastran.bdf.bdf_interface.pybdf import (
    BDFInputPy, _show_bad_file)


def read_bdf(bdf_filename=None, validate=True, xref=True, punch=False,
             encoding=None, log=None, debug=True, mode='msc'):
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
    validate : bool
        runs various checks on the BDF (default=True)
    xref :  bool
        should the bdf be cross referenced (default=True)
    punch : bool
        indicates whether the file is a punch file (default=False)
    encoding : str
        the unicode encoding (default=None; system default)

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
    model.read_bdf(bdf_filename=bdf_filename, validate=validate,
                   xref=xref, punch=punch, read_includes=True, encoding=encoding)

    #if 0:
        #keys_to_suppress = []
        #method_names = model.object_methods(keys_to_skip=keys_to_suppress)

        #methods_to_remove = [
            #'_process_card', 'read_bdf', 'disable_cards', 'set_dynamic_syntax',
            #'create_card_object', 'create_card_object_fields', 'create_card_object_list',

            #'add_AECOMP', 'add_AEFACT', 'add_AELINK', 'add_AELIST', 'add_AEPARM', 'add_AERO',
            #'add_AEROS', 'add_AESTAT', 'add_AESURF', 'add_ASET', 'add_BCRPARA', 'add_BCTADD',
            #'add_BCTPARA', 'add_BCTSET', 'add_BSET', 'add_BSURF', 'add_BSURFS', 'add_CAERO',
            #'add_DIVERG',
            #'add_CSET', 'add_CSSCHD', 'add_DAREA', 'add_DCONADD', 'add_DCONSTR', 'add_DDVAL',
            #'add_DELAY', 'add_DEQATN', 'add_DESVAR', 'add_DLINK', 'add_DMI', 'add_DMIG',
            #'add_DMIJ', 'add_DMIJI', 'add_DMIK', 'add_DPHASE', 'add_DRESP', 'add_DTABLE',
            #'add_DVMREL', 'add_DVPREL', 'add_EPOINT', 'add_FLFACT', 'add_FLUTTER', 'add_FREQ',
            #'add_GUST', 'add_LSEQ', 'add_MKAERO', 'add_MONPNT', 'add_NLPARM', 'add_NLPCI',
            #'add_PAERO', 'add_PARAM', 'add_PBUSHT', 'add_PDAMPT', 'add_PELAST', 'add_PHBDY',
            #'add_QSET', 'add_SEBSET', 'add_SECSET', 'add_SEQSET', 'add_SESET', 'add_SET',
            #'add_SEUSET', 'add_SPLINE', 'add_spoint', 'add_tempd', 'add_TF', 'add_TRIM',
            #'add_TSTEP', 'add_TSTEPNL', 'add_USET',

            #'add_card', 'add_card_fields', 'add_cmethod', 'add_constraint',
            #'add_constraint_MPC', 'add_constraint_MPCADD',
            #'add_constraint_SPC', 'add_constraint_SPCADD',
            #'add_convection_property', 'add_coord', 'add_creep_material', 'add_damper',
            #'add_dload', '_add_dload_entry', 'add_element', 'add_hyperelastic_material',
            #'add_load', 'add_mass', 'add_material_dependence', 'add_method', 'add_node',
            #'add_plotel', 'add_property', 'add_property_mass', 'add_random_table',
            #'add_rigid_element', 'add_structural_material', 'add_suport', 'add_suport1',
            #'add_table', 'add_table_sdamping', 'add_thermal_BC', 'add_thermal_element',
            #'add_thermal_load', 'add_thermal_material',

            #'set_as_msc',
            #'set_as_nx',

            #'pop_parse_errors',
            ##'pop_xref_errors',
            #'set_error_storage',
            #'is_reject',
        #]
        #for method_name in method_names:
            #if method_name not in methods_to_remove + keys_to_suppress:
                ##print(method_name)
                #pass
            #else:
                ### TODO: doesn't work...
                ##delattr(model, method_name)
                #pass
        #model.get_bdf_stats()
    return model


class BDF(AddCard, CrossReference, WriteMesh, GetMethods):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    #__slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug: str | bool | None,
                 log: Optional[SimpleLogger]=None, mode: str='msc'):
        """
        Initializes the BDF object

        Parameters
        ----------
        debug : bool/None
            used to set the logger if no logger is passed in
                True:  logs debug/info/error messages
                False: logs info/error messages
                None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        """
        AddCard.__init__(self)
        CrossReference.__init__(self)
        WriteMesh.__init__(self)
        GetMethods.__init__(self)
        assert debug in [True, False, None], 'debug=%r' % debug
        self.echo = False
        self.read_includes = True

        # file management parameters
        self.active_filenames = []
        self.active_filename = None
        self.include_dir = ''
        self.dumplines = False

        # this flag will be flipped to True someday (and then removed), but
        # doesn't support 100% of cards yet.  It enables a new method for card
        # parsing.
        #
        # 80.3 seconds -> 67.2 seconds for full_bay model
        # (multiple BDF passes among other things)
        self._fast_add = True

        self.log = get_logger(log=log, level=debug, **log_args)

        #: list of all read in cards - useful in determining if entire BDF
        #: was read & really useful in debugging
        self.card_count = {}
        #: stores the card_count of cards that have been rejected
        self.reject_count = {}

        #: was an ENDDATA card found
        #self.foundEndData = False

        #: useful in debugging errors in input
        self.debug = debug

        #: flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False

        #: lines that were rejected b/c they were for a card that isn't supported
        self.reject_lines = []

        #: cards that were created, but not processed
        self.reject_cards = []

        # self.__init_attributes()

        #: the list of possible cards that will be parsed
        self.cards_to_read = set([
            'GRID', 'SPOINT', 'EPOINT', 'POINT', 'POINTAX',
            'PARAM', ## params

            # coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',

            'PELAS', 'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',

            'CROD', 'PROD', 'CONROD',
            'CTUBE', 'PTUBE',

            'PBAR', 'PBARL', 'CBAR',
            'CBEAM',

            'PSHEAR', 'CSHEAR',

            'CQUAD4', 'CTRIA3', 'CQUAD8', 'CTRIA6',
            'PSHELL', 'PCOMP', 'PCOMPG',

            'PSOLID', 'PLSOLID',
            'CTETRA', 'CTETRA4', 'CTETRA10',
            'CPYRAM', 'CPYRAM5', 'CPYRAM13',
            'CPENTA', 'CPENTA6', 'CPENTA15',
            'CHEXA', 'CHEXA8', 'CHEXA20',

            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            #'PBUSH', 'PBUSH1D', 'PBUSH2D',

            'CONM1', 'CONM2',
            'PLOTEL',
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE',

            'MAT1', 'MAT8',

            # loads
            'LOAD', 'GRAV',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'PLOAD', 'PLOAD2', 'PLOAD4', 'PLOADX1',
            'TLOAD1', 'TLOAD2', 'DELAY',
            'RLOAD1', 'DPHASE', #'RLOAD2',


            # constraints
            'SPC', 'SPCADD', 'SPC1', 'SPCD',
            'MPC', 'MPCADD',

            # aero cards
            'AERO',  ## aero
            'AEROS',  ## aeros
            'GUST',  ## gusts
            'FLUTTER',   ## flutters
            'FLFACT',   ## flfacts
            'MKAERO1', 'MKAERO2',  ## mkaeros
            'AECOMP',   ## aecomps
            'AEFACT',   ## aefacts
            'AELINK',   ## aelinks
            'AELIST',   ## aelists
            'AEPARM',   ## aeparams
            'AESTAT',   ## aestats
            'AESURF',  ## aesurf
            #'AESURFS', ## aesurfs
            'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4',  ## caeros
            # 'CAERO5',
            'PAERO1', 'PAERO2', 'PAERO3',  ## paeros
            'PAERO4', # 'PAERO5',
            'MONPNT1',                                   ## monitor_points
            'SPLINE1', 'SPLINE2', 'SPLINE4', 'SPLINE5',  ## splines
            #'SPLINE3', 'SPLINE6', 'SPLINE7',
            'TRIM',  ## trims
            'CSSCHD', ## csschds
            'DIVERG', ## divergs

            # ---- dynamic cards ---- #
            'DAREA',  ## dareas
            'DPHASE',  ## dphases
            'DELAY',  ## delays
            'NLPARM',  ## nlparms
            'ROTORG', 'ROTORD', ## rotors
            'NLPCI',  ## nlpcis
            'TSTEP',  ## tsteps
            'TSTEPNL', 'TSTEP1',  ## tstepnls
            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',

            # optimization cards
            'DEQATN', 'DTABLE',
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
            'DVCREL1', 'DVCREL2',
            'DVPREL1', 'DVPREL2',
            'DVMREL1', 'DVMREL2',
            'DOPTPRM', 'DLINK', 'DCONADD', 'DVGRID',

            'SET1', 'SET3',  ## sets
            'ASET', 'ASET1',  ## asets
            'BSET', 'BSET1',  ## bsets
            'CSET', 'CSET1',  ## csets
            'QSET', 'QSET1',  ## qsets
            'USET', 'USET1',  ## usets

            ## suport/suport1/se_suport
            'SUPORT', 'SUPORT1', 'SESUP',

            #: methods
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods
            'EIGC', 'EIGP',

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ])

        case_control_cards = {'FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                              'TEMP', 'TSTEPNL', 'INCLUDE'}
        self._unique_bulk_data_cards = self.cards_to_read.difference(case_control_cards)

        #: / is the delete from restart card
        self.special_cards = ['DEQATN', '/']
        self._make_card_parser()
        if self.is_msc:
            self.set_as_msc()
        elif self.is_nx:
            self.set_as_nx()
        #elif self.is_optistruct:
            #self.set_as_optistruct()
        #elif self.is_radioss:
            #self.set_as_radioss()
        else:
            msg = 'mode=%r is not supported; modes=[msc, nx]' % self._nastran_format
            raise NotImplementedError(msg)

    def __getstate__(self):
        """clears out a few variables in order to pickle the object"""
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        # Remove the unpicklable entries.
        #del state['spcObject'], state['mpcObject'],
        del state['_card_parser'], state['log']
        return state

    def save(self, obj_filename='model.obj', unxref=True):
        """
        ..warning:: doesn't work right
        """
        #del self.log
        #del self.spcObject
        #del self.mpcObject
        #del self._card_parser, self._card_parser_prepare

        #try:
            #del self.log
        #except AttributeError:
            #pass
        #self.case_control_lines = str(self.case_control_deck).split('\n')
        #del self.case_control_deck

        if unxref:
            self.uncross_reference()
        with open(obj_filename, 'w') as obj_file:
            dump(self, obj_file)

    def load(self, obj_filename='model.obj'):
        """
        ..warning:: doesn't work right
        """
        return
        #del self.log
        #del self.spcObject
        #del self.mpcObject
        #lines = print(self.case_control_deck)
        #self.case_control_lines = lines.split('\n')
        #del self.case_control_deck
        #self.uncross_reference()
        #import types
        with open(obj_filename, "r") as obj_file:
            obj = load(obj_file)

        keys_to_skip = [
            'case_control_deck',
            'log', #'mpcObject', 'spcObject',
            'node_ids', 'coord_ids', 'element_ids', 'property_ids',
            'material_ids', 'caero_ids', 'is_long_ids',
            'nnodes', 'ncoords', 'nelements', 'nproperties',
            'nmaterials', 'ncaeros',

            'point_ids', 'subcases',
            '_card_parser',
        ]
        for key in object_attributes(self, mode="all", keys_to_skip=keys_to_skip):
            if key.startswith('__') and key.endswith('__'):
                continue

            #print('key =', key)
            val = getattr(obj, key)
            #print(key)
            #if isinstance(val, types.FunctionType):
                #continue
            setattr(self, key, val)

        self.case_control_deck = CaseControlDeck(self.case_control_lines, log=self.log)
        self.log.debug('done loading!')

    def replace_cards(self, replace_model):
        """
        Replaces the common cards from the current (self) model from the
        ones in the new replace_model.  The intention is that you're
        going to replace things like PSHELLs and DESVARs from a pch file
        in order to update your BDF with the optimized geometry.

        .. todo:: only does a subset of cards.

        Notes
        -----
        loads/spcs (not supported) are tricky because you
        can't replace cards one-to-one...not sure what to do.
        """
        for nid, node in replace_model.nodes.items():
            self.nodes[nid] = node
        for eid, elem in replace_model.elements.items():
            self.elements[eid] = elem
        for eid, elem in replace_model.rigid_elements.items():
            self.rigid_elements[eid] = elem
        for pid, prop in replace_model.properties.items():
            self.properties[pid] = prop
        for mid, mat in replace_model.materials.items():
            self.materials[mid] = mat

        for dvid, desvar in replace_model.desvars.items():
            self.desvars[dvid] = desvar
        for dvid, dvprel in replace_model.dvprels.items():
            self.dvprels[dvid] = dvprel
        for dvid, dvmrel in replace_model.dvmrels.items():
            self.dvmrels[dvid] = dvmrel
        for dvid, dvgrid in replace_model.dvgrids.items():
            self.dvgrids[dvid] = dvgrid

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        Parameters
        ----------
        cards : list[str]; set[str]
            a list/set of cards that should not be read

        .. python ::

            bdfModel.disable_cards(['DMIG', 'PCOMP'])
        """
        if cards is None:
            return
        elif isinstance(cards, str):
            disable_set = set([cards])
        else:
            disable_set = set(cards)
        self.cards_to_read = self.cards_to_read.difference(disable_set)

    def set_error_storage(self, nparse_errors: int=100,
                          stop_on_parsing_error: bool=True,
                          nxref_errors: int=100,
                          stop_on_xref_error: bool=True):
        """
        Catch parsing errors and store them up to print them out all at once
        (not all errors are caught).

        Parameters
        ----------
        nparse_errors : int
            how many parse errors should be stored
            (default=0; all=None; no storage=0)
        stop_on_parsing_error : bool
            should an error be raised if there
            are parsing errors (default=True)
        nxref_errors : int
            how many cross-reference errors
            should be stored (default=0; all=None; no storage=0)
        stop_on_xref_error : bool
            should an error be raised if there
            are cross-reference errors (default=True)
        """
        self.xref_obj.set_error_storage(
            nparse_errors=nparse_errors,
            stop_on_parsing_error=stop_on_parsing_error,
            nxref_errors=nxref_errors,
            stop_on_xref_error=stop_on_xref_error,
        )
        assert isinstance(nparse_errors, int), type(nparse_errors)
        #assert isinstance(nxref_errors, int), type(nxref_errors)
        self._nparse_errors = nparse_errors
        #self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        #self._stop_on_xref_error = stop_on_xref_error

    def validate(self):
        """runs some checks on the input data beyond just type checking"""
        return
        for nid, node in sorted(self.nodes.items()):
            node.validate()
        for cid, coord in sorted(self.coords.items()):
            coord.validate()
        for eid, elem in sorted(self.elements.items()):
            elem.validate()
        for pid, prop in sorted(self.properties.items()):
            prop.validate()

        for eid, elem in sorted(self.rigid_elements.items()):
            elem.validate()
        for eid, plotel in sorted(self.plotels.items()):
            plotel.validate()
        #for eid, mass in sorted(self.masses.items()):
            #mass.validate()
        for pid, property_mass in sorted(self.properties_mass.items()):
            property_mass.validate()

        #------------------------------------------------
        for mid, mat in sorted(self.materials.items()):
            mat.validate()
        for mid, mat in sorted(self.thermal_materials.items()):
            mat.validate()
        for mid, mat in sorted(self.MATS1.items()):
            mat.validate()
        for mid, mat in sorted(self.MATS3.items()):
            mat.validate()
        for mid, mat in sorted(self.MATS8.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT1.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT2.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT3.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT4.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT5.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT8.items()):
            mat.validate()
        for mid, mat in sorted(self.MATT9.items()):
            mat.validate()
        for mid, mat in sorted(self.creep_materials.items()):
            mat.validate()
        for mid, mat in sorted(self.hyperelastic_materials.items()):
            mat.validate()

        #------------------------------------------------
        for key, loads in sorted(self.loads.items()):
            for loadi in loads:
                loadi.validate()
        for key, tic in sorted(self.tics.items()):
            tic.validate()
        for key, dloads in sorted(self.dloads.items()):
            for dload in dloads:
                dload.validate()
        for key, dload_entries in sorted(self.dload_entries.items()):
            for dload_entry in dload_entries:
                dload_entry.validate()

        #------------------------------------------------
        for key, nlpci in sorted(self.nlpcis.items()):
            nlpci.validate()
        for key, nlparm in sorted(self.nlparms.items()):
            nlparm.validate()
        for key, tstep in sorted(self.tsteps.items()):
            tstep.validate()
        for key, tstepnl in sorted(self.tstepnls.items()):
            tstepnl.validate()
        for key, transfer_functions in sorted(self.transfer_functions.items()):
            for transfer_function in transfer_functions:
                transfer_function.validate()
        for key, delay in sorted(self.delays.items()):
            delay.validate()

        #------------------------------------------------
        if self.aeros is not None:
            self.aeros.validate()
        for caero_id, caero in sorted(self.caeros.items()):
            caero.validate()
        for key, paero in sorted(self.paeros.items()):
            paero.validate()
        for spline_id, spline in sorted(self.splines.items()):
            spline.validate()

        for key, aecomp in sorted(self.aecomps.items()):
            aecomp.validate()
        for key, aefact in sorted(self.aefacts.items()):
            aefact.validate()
        for key, aelinks in sorted(self.aelinks.items()):
            for aelink in aelinks:
                aelink.validate()
        for key, aeparam in sorted(self.aeparams.items()):
            aeparam.validate()
        for key, aesurf in sorted(self.aesurf.items()):
            aesurf.validate()
        for key, aesurfs in sorted(self.aesurfs.items()):
            aesurfs.validate()
        for key, aestat in sorted(self.aestats.items()):
            aestat.validate()
        for key, trim in sorted(self.trims.items()):
            trim.validate()
        for key, diverg in sorted(self.divergs.items()):
            diverg.validate()
        for key, csschd in sorted(self.csschds.items()):
            csschd.validate()
        for monitor in self.monitor_points:
            monitor.validate()

        #------------------------------------------------
        if self.aero is not None:
            self.aero.validate()
        for key, flfact in sorted(self.flfacts.items()):
            flfact.validate()
        for key, flutter in sorted(self.flutters.items()):
            flutter.validate()
        for key, gust in sorted(self.gusts.items()):
            gust.validate()
        #self.mkaeros = []

        #------------------------------------------------
        for key, bcs in sorted(self.bcs.items()):
            for bc in bcs:
                bc.validate()
        for key, phbdy in sorted(self.phbdys.items()):
            phbdy.validate()
        for key, convection_property in sorted(self.convection_properties.items()):
            convection_property.validate()
        for key, tempd in sorted(self.tempds.items()):
            tempd.validate()
        #------------------------------------------------
        for key, bcrpara in sorted(self.bcrparas.items()):
            bcrpara.validate()
        for key, bctadd in sorted(self.bctadds.items()):
            bctadd.validate()
        for key, bctpara in sorted(self.bctparas.items()):
            bctpara.validate()
        for key, bctset in sorted(self.bctsets.items()):
            bctset.validate()
        for key, bsurf in sorted(self.bsurf.items()):
            bsurf.validate()
        for key, bsurfs in sorted(self.bsurfs.items()):
            bsurfs.validate()

        #------------------------------------------------
        for key, suport1 in sorted(self.suport1.items()):
            suport1.validate()
        for suport in self.suport:
            suport.validate()
        for se_suport in self.se_suport:
            se_suport.validate()

        for key, spcs in sorted(self.spcs.items()):
            for spc in spcs:
                spc.validate()
        for key, spcadd in sorted(self.spcadds.items()):
            spcadd.validate()

        for key, mpcs in sorted(self.mpcs.items()):
            for mpc in mpcs:
                mpc.validate()
        for key, mpcadd in sorted(self.mpcadds.items()):
            mpcadd.validate()

        #------------------------------------------------
        #for key, darea in sorted(self.dareas.items()):
            #darea.validate()
        #for key, dphase in sorted(self.dphases.items()):
            #dphase.validate()

        for pid, pbusht in sorted(self.pbusht.items()):
            pbusht.validate()
        for pid, pdampt in sorted(self.pdampt.items()):
            pdampt.validate()
        for pid, pelast in sorted(self.pelast.items()):
            pelast.validate()

        for pid, frequency in sorted(self.frequencies.items()):
            frequency.validate()
        #------------------------------------------------
        for key, dmi in sorted(self.dmi.items()):
            dmi.validate()
        for key, dmig in sorted(self.dmig.items()):
            dmig.validate()
        for key, dmij in sorted(self.dmij.items()):
            dmij.validate()
        for key, dmiji in sorted(self.dmiji.items()):
            dmiji.validate()
        for key, dmik in sorted(self.dmik.items()):
            dmik.validate()
        #------------------------------------------------
        #self.asets = []
        #self.bsets = []
        #self.csets = []
        #self.qsets = []
        #self.usets = {}

        ##: SExSETy
        #self.se_bsets = []
        #self.se_csets = []
        #self.se_qsets = []
        #self.se_usets = {}
        #self.se_sets = {}

        for key, sets in sorted(self.sets.items()):
            sets.validate()
        for key, uset in sorted(self.usets.items()):
            for useti in uset:
                useti.validate()

        for aset in self.asets:
            aset.validate()
        for omit in self.omits:
            omit.validate()
        for bset in self.bsets:
            bset.validate()
        for cset in self.csets:
            cset.validate()
        for qset in self.qsets:
            qset.validate()

        for key, se_set in sorted(self.se_sets.items()):
            se_set.validate()
        for key, se_uset in sorted(self.se_usets.items()):
            se_uset.validate()
        for se_bset in self.se_bsets:
            se_bset.validate()
        for se_cset in self.se_csets:
            se_cset.validate()
        for se_qset in self.se_qsets:
            se_qset.validate()
        #------------------------------------------------
        for key, table in sorted(self.tables.items()):
            table.validate()
        for key, table in sorted(self.tables_d.items()):
            table.validate()
        for key, table in sorted(self.tables_m.items()):
            table.validate()
        for key, random_table in sorted(self.random_tables.items()):
            random_table.validate()
        for key, table_sdamping in sorted(self.tables_sdamping.items()):
            table_sdamping.validate()
        #------------------------------------------------
        for key, method in sorted(self.methods.items()):
            method.validate()
        for key, cmethod in sorted(self.cMethods.items()):
            cmethod.validate()
        #------------------------------------------------
        for key, dconadd in sorted(self.dconadds.items()):
            dconadd.validate()
        for key, dconstrs in sorted(self.dconstrs.items()):
            for dconstr in dconstrs:
                dconstr.validate()
        for key, desvar in sorted(self.desvars.items()):
            desvar.validate()
        for key, ddval in sorted(self.ddvals.items()):
            ddval.validate()
        for key, dlink in sorted(self.dlinks.items()):
            dlink.validate()
        for key, dresp in sorted(self.dresps.items()):
            dresp.validate()

        if self.dtable is not None:
            self.dtable.validate()
        if self.doptprm is not None:
            self.doptprm.validate()
        for key, dequation in sorted(self.dequations.items()):
            dequation.validate()
        for key, dvprel in sorted(self.dvprels.items()):
            dvprel.validate()
        for key, dvmrel in sorted(self.dvmrels.items()):
            dvmrel.validate()
        for key, dvcrel in sorted(self.dvcrels.items()):
            dvcrel.validate()
        for key, dscreen in sorted(self.dscreen.items()):
            dscreen.validate()
        for dvid, dvgrid in self.dvgrids.items():
            dvgrid.validate()
        #------------------------------------------------

    def read_bdf(self, bdf_filename=None,
                 validate=True, xref=True, punch=False, read_includes=True, encoding=None):
        """
        Read method for the bdf files

        Parameters
        ----------
        bdf_filename : str / None
            the input bdf (default=None; popup a dialog)
        validate : bool
            runs various checks on the BDF (default=True)
        xref :  bool
            should the bdf be cross referenced (default=True)
        punch : bool
            indicates whether the file is a punch file (default=False)
        read_includes : bool
            indicates whether INCLUDE files should be read (default=True)
        encoding : str
            the unicode encoding (default=None; system default)

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

        """
        self._read_bdf_helper(bdf_filename, encoding, punch, read_includes)

        self._parse_primary_file_header(bdf_filename)

        self.log.debug(f'---starting BDF.read_bdf of {self.bdf_filename}---')

        #executive_control_lines, case_control_lines, \
            #bulk_data_lines = self.get_lines(self.bdf_filename, self.punch)
        obj = BDFInputPy(self.read_includes, self.dumplines, self._encoding,
                         log=self.log, debug=self.debug)
        out = obj.get_lines(bdf_filename, punch=self.punch)
        #system_lines, executive_control_lines, case_control_lines, bulk_data_lines = out
        (system_lines, executive_control_lines, case_control_lines,
         bulk_data_lines, bulk_data_ilines,
         additional_deck_lines, additional_deck_ilines) = out
        self._set_pybdf_attributes(obj)

        self.case_control_lines = case_control_lines
        self.executive_control_lines = executive_control_lines

        sol, method, sol_iline, app = parse_executive_control_deck(executive_control_lines)
        self.update_solution(sol, method, sol_iline)

        self.case_control_deck = CaseControlDeck(self.case_control_lines, self.log)
        #print(self.object_attributes())
        self.case_control_deck.solmap_to_value = self._solmap_to_value
        self.case_control_deck.rsolmap_to_str = self.rsolmap_to_str

        if self._is_cards_dict:
            cards, card_count = self.get_bdf_cards_dict(bulk_data_lines)
        else:
            cards, card_count = self.get_bdf_cards(bulk_data_lines)
        self._parse_cards(cards, card_count)

        if 0 and self.values_to_skip:
            for key, values in self.values_to_skip.items():
                dict_values = getattr(self, key)
                if not isinstance(dict_values, dict):
                    msg = '%r is an invalid type; only dictionaries are supported' % key
                    raise TypeError(msg)
                for value in values:
                    del dict_values[value]
            # TODO: redo get_card_ids_by_card_types & card_count

        #self.pop_parse_errors()
        fill_dmigs(self)

        if validate:
            self.validate()

        self.cross_reference(xref=xref)
        self._xref = xref

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
        #self.pop_xref_errors()

    def _set_pybdf_attributes(self, obj):
        """common method for all functions that use BDFInputPy"""
        #self.reject_lines += obj.reject_lines
        self.active_filenames += obj.active_filenames
        self.active_filename = obj.active_filename
        self.include_dir = obj.include_dir

    def _read_bdf_helper(self, bdf_filename: Optional[PathLike],
                         encoding: Optional[str],
                         punch: bool, read_includes: bool):
        """creates the file loading if bdf_filename is None"""
        #self.set_error_storage(nparse_errors=None, stop_on_parsing_error=True,
        #                       nxref_errors=None, stop_on_xref_error=True)
        if encoding is None:
            encoding = sys.getdefaultencoding()
        self._encoding = encoding
        if bdf_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran BDF (*.bdf; *.dat; *.nas; *.pch, *.ecd)|" \
                "*.bdf;*.dat;*.nas;*.pch|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran BDF (*.bdf *.dat *.nas *.pch *.ecd);;All files (*)"
            title = 'Please select a BDF/DAT/PCH/ECD to load'
            bdf_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)[0]
            assert bdf_filename is not None, bdf_filename

        check_path(bdf_filename, 'bdf_filename')
        if bdf_filename.lower().endswith('.pch'):  # .. todo:: should this be removed???
            punch = True

        #: the active filename (string)
        self.bdf_filename = bdf_filename

        #: is this a punch file (no executive control deck)
        self.punch = punch
        self.read_includes = read_includes
        self.active_filenames = []

    def fill_dmigs(self):
        """fills the DMIx cards with the column data that's been stored"""
        return
        #for name, card_comments in self._dmig_temp.items():
            #card0, comment0 = card_comments[0]
            #card_name = card0[0]
            #card_name = card_name.rstrip(' *').upper()

            #if card_name == 'DMIG':
                ## if field2 == 'UACCEL':  # special DMIG card
                #card = self.dmig[name]
            #elif card_name == 'DMI':
                #card = self.dmi[name]
            #elif card_name == 'DMIJ':
                #card = self.dmij[name]
            #elif card_name == 'DMIJI':
                #card = self.dmiji[name]
            #elif card_name == 'DMIK':
                #card = self.dmik[name]
            #else:
                #raise NotImplementedError(card_name)

            #for (card_obj, comment) in card_comments:
                #card._add_column(card_obj, comment=comment)
            #card.finalize()

        #self._dmig_temp = defaultdict(list)

    def pop_parse_errors(self):
        """raises an error if there are parsing errors"""
        xref_obj = self.xref_obj
        if self._stop_on_parsing_error:
            if self._iparse_errors == 1 and self._nparse_errors == 0:
                raise
            is_error = False
            msg = ''
            if self._duplicate_elements:
                duplicate_eids = [elem.eid for elem in self._duplicate_elements]
                uduplicate_eids = np.unique(duplicate_eids)
                msg += 'self.elements IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_element=\n%s\n' % str(self.elements[eid])
                    msg += 'new_elements=\n'
                    for elem, eidi in zip(self._duplicate_elements, duplicate_eids):
                        if eidi == eid:
                            msg += str(elem)
                    msg += '\n'
                    is_error = True
                    raise DuplicateIDsError(msg)

            if self._duplicate_properties:
                duplicate_pids = [prop.pid for prop in self._duplicate_properties]
                uduplicate_pids = np.unique(duplicate_pids)
                msg += 'self.properties IDs are not unique=%s\n' % uduplicate_pids
                for pid in duplicate_pids:
                    msg += 'old_property=\n%s\n' % str(self.properties[pid])
                    msg += 'new_properties=\n'
                    for prop, pidi in zip(self._duplicate_properties, duplicate_pids):
                        if pidi == pid:
                            msg += str(prop)
                    msg += '\n'
                    is_error = True

            if self._duplicate_masses:
                duplicate_eids = [elem.eid for elem in self._duplicate_masses]
                uduplicate_eids = np.unique(duplicate_eids)
                msg += 'self.massses IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_mass=\n%s\n' % str(self.masses[eid])
                    msg += 'new_masses=\n'
                    for elem, eidi in zip(self._duplicate_masses, duplicate_eids):
                        if eidi == eid:
                            msg += str(elem)
                    msg += '\n'
                    is_error = True

            if self._duplicate_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_material=\n%s\n' % str(self.materials[mid])
                    msg += 'new_materials=\n'
                    for mat, midi in zip(self._duplicate_materials, duplicate_mids):
                        if midi == mid:
                            msg += str(mat)
                    msg += '\n'
                    is_error = True

            if self._duplicate_thermal_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_thermal_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.thermal_materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_thermal_material=\n%s\n' % str(self.thermal_materials[mid])
                    msg += 'new_thermal_materials=\n'
                    for mat, midi in zip(self._duplicate_thermal_materials, duplicate_mids):
                        if midi == mid:
                            msg += str(mat)
                    msg += '\n'
                    is_error = True

            if self._duplicate_coords:
                duplicate_cids = [coord.cid for coord in self._duplicate_coords]
                uduplicate_cids = np.unique(duplicate_cids)
                msg += 'self.coords IDs are not unique=%s\n' % uduplicate_cids
                for cid in uduplicate_cids:
                    msg += 'old_coord=\n%s\n' % str(self.coords[cid])
                    msg += 'new_coords=\n'
                    for coord, cidi in zip(self._duplicate_coords, duplicate_cids):
                        if cidi == cid:
                            msg += str(coord)
                    msg += '\n'
                    is_error = True

            if is_error:
                msg = 'There are duplicate cards.\n\n' + msg

            if xref_obj._stop_on_xref_error:
                msg += 'There are parsing errors.\n\n'
                for (card, an_error) in self._stored_parse_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    msg += 'xref error: %s\n\n'% an_error[0]
                    is_error = True

            if is_error:
                self.log.error('%s' % msg)
                raise DuplicateIDsError(msg.rstrip())

    def pop_xref_errors(self):
        """raises an error if there are cross-reference errors"""
        self.xref_obj.pop_xref_errors()

    def get_bdf_cards(self, bulk_data_lines):
        """Parses the BDF lines into a list of card_lines"""
        cards = []
        #cards = defaultdict(list)
        card_count = defaultdict(int)
        full_comment = ''
        card_lines = []
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)
        for i, line in enumerate(bulk_data_lines):
            #print('    backup=%r' % backup_comment)
            comment = ''
            if '$' in line:
                line, comment = line.split('$', 1)
            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    if self.echo:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    # cards[old_card_name].append([full_comment, card_lines])

                    # new list version
                    #if full_comment:
                        #print('full_comment = ', full_comment)
                    cards.append([old_card_name, _prep_comment(full_comment), card_lines])

                    card_count[old_card_name] += 1
                    card_lines = []
                    full_comment = ''

                    if old_card_name == 'ECHOON':
                        self.echo = True
                    elif old_card_name == 'ECHOOFF':
                        self.echo = False
                old_card_name = card_name.rstrip(' *')
                if old_card_name == 'ENDDATA':
                    self.card_count['ENDDATA'] = 1
                    if nlines - i > 1:
                        nleftover = nlines - i - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)
            if line.rstrip():
                card_lines.append(line)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += comment + '\n'
                #print('add backup=%r' % backup_comment)
            #elif comment:
                #backup_comment += '$' + comment + '\n'

        if card_lines:
            if self.echo:
                self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)

            # old dictionary version
            #cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            #if backup_comment + full_comment:
                #print('backup_comment + full_comment = ', backup_comment + full_comment)
            cards.append([old_card_name, _prep_comment(backup_comment + full_comment), card_lines])
            card_count[old_card_name] += 1
        return cards, card_count

    def get_bdf_cards_dict(self, bulk_data_lines):
        """Parses the BDF lines into a list of card_lines"""
        cards = defaultdict(list)
        card_count = defaultdict(int)
        full_comment = ''
        card_lines = []
        old_card_name = None
        backup_comment = ''
        nlines = len(bulk_data_lines)
        for i, line in enumerate(bulk_data_lines):
            #print('    backup=%r' % backup_comment)
            comment = ''
            if '$' in line:
                line, comment = line.split('$', 1)
            card_name = line.split(',', 1)[0].split('\t', 1)[0][:8].rstrip().upper()
            if card_name and card_name[0] not in ['+', '*']:
                if old_card_name:
                    if self.echo:
                        self.log.info('Reading %s:\n' %
                                      old_card_name + full_comment + ''.join(card_lines))

                    # old dictionary version
                    cards[old_card_name].append([full_comment, card_lines])

                    # new list version
                    #cards.append([old_card_name, full_comment, card_lines])

                    card_count[old_card_name] += 1
                    card_lines = []
                    full_comment = ''

                    if old_card_name == 'ECHOON':
                        self.echo = True
                    elif old_card_name == 'ECHOOFF':
                        self.echo = False
                old_card_name = card_name.rstrip(' *')
                if old_card_name == 'ENDDATA':
                    self.card_count['ENDDATA'] = 1
                    if nlines - i > 1:
                        nleftover = nlines - i - 1
                        msg = 'exiting due to ENDDATA found with %i lines left' % nleftover
                        self.log.debug(msg)
                    return cards, card_count
                #print("card_name = %s" % card_name)

            comment = _clean_comment(comment)
            if line.rstrip():
                card_lines.append(line)
                if backup_comment:
                    if comment:
                        full_comment += backup_comment + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += comment + '\n'
                #print('add backup=%r' % backup_comment)
            #elif comment:
                #backup_comment += comment + '\n'

        if card_lines:
            if self.echo:
                self.log.info('Reading %s:\n' % old_card_name + full_comment + ''.join(card_lines))
            #print('end_add %s' % card_lines)

            # old dictionary version
            cards[old_card_name].append([backup_comment + full_comment, card_lines])

            # new list version
            #cards.append([old_card_name, backup_comment + full_comment, card_lines])
            card_count[old_card_name] += 1
        return cards, card_count

    def update_solution(self, sol, method, sol_iline):
        """
        Updates the overall solution type (e.g. 101,200,600)

        Parameters
        ----------
        sol : int
            the solution type (101, 103, etc)
        method : str
            the solution method (only for SOL=600)
        sol_iline : int
            the line to put the SOL/method on
        """
        self.sol_iline = sol_iline
        # the integer of the solution type (e.g. SOL 101)
        if sol is None:
            self.sol = None
            self.sol_method = None
            return

        try:
            self.sol = int(sol)
        except ValueError:
            try:
                self.sol = self._solmap_to_value[sol]
            except KeyError:
                self.sol = sol

        if self.sol == 600:
            #: solution 600 method modifier
            self.sol_method = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.sol_method))
        else:  # very common
            self.sol_method = None

    def set_dynamic_syntax(self, dict_of_vars):
        """
        Uses the OpenMDAO syntax of %varName in an embedded BDF to
        update the values for an optimization study.

        Parameters
        ----------
        dict_of_vars : dict[str] = int/float/str
            dictionary of 7 character variable names to map.

        .. code-block:: python

           GRID, 1, %xVar, %yVar, %zVar

           >>> dict_of_vars = {'xVar': 1.0, 'yVar', 2.0, 'zVar':3.0}
           >>> bdf = BDF()
           >>> bdf.set_dynamic_syntax(dict_of_vars)
           >>> bdf,read_bdf(bdf_filename, xref=True)

        Notes
        -----
        Case sensitivity is supported.

        Variables should be 7 characters or less to fit in an
        8-character field.

        .. warning:: Type matters!
        """
        self.dict_of_vars = {}
        assert len(dict_of_vars) > 0, 'nvars = %s' % len(dict_of_vars)
        for (key, value) in sorted(dict_of_vars.items()):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) >= 1, ('min length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            if not isinstance(key, str):
                msg = 'key=%r must be a string.  type=%s' % (key, type(key))
                raise TypeError(msg)
            self.dict_of_vars[key] = value
        self._is_dynamic_syntax = True

    def is_reject(self, card_name):
        """
        Can the card be read.

        If the card is rejected, it's added to self.reject_count

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'
        """
        if card_name.startswith('='):
            return False
        elif card_name in self.cards_to_read:
            return False
        if card_name:
            if card_name not in self.reject_count:
                self.reject_count[card_name] = 0
            self.reject_count[card_name] += 1
        return True

    def _process_card(self, card_lines):
        """
        Converts card_lines into a card.
        Considers dynamic syntax and removes empty fields

        Parameters
        ----------
        card_lines : list[str]
            list of strings that represent the card's lines

        Returns
        -------
        fields : list[str]
            the parsed card's fields
        card_name : str
            the card's name

        .. code-block:: python

           >>> card_lines = ['GRID,1,,1.0,2.0,3.0,,']
           >>> model = BDF()
           >>> fields, card_name = model._process_card(card_lines)
           >>> fields
           ['GRID', '1', '', '1.0', '2.0', '3.0']
           >>> card_name
           'GRID'

        """
        card_name = self._get_card_name(card_lines)
        fields = to_fields(card_lines, card_name)
        if self._is_dynamic_syntax:
            fields = [self._parse_dynamic_syntax(field) if '%' in
                      field[0:1] else field for field in fields]
        card = wipe_empty_fields(fields)
        card[0] = card_name
        return card

    def create_card_object(self, card_lines, card_name, is_list=True, has_none=True):
        """
        Creates a BDFCard object, which is really just a list that
        allows indexing past the last field

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is list of card_lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']
        card_name : str
            the card_name -> 'GRID'
        is_list : bool; default=True
            True : this is a list of fields
            False : this is a list of lines
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card
        card : list[str]
            the card with empty fields removed
        """
        card_name = card_name.upper()
        self.increase_card_count(card_name)
        if card_name in ['DEQATN', 'PBRSECT', 'PBMSECT']:
            card_obj = card_lines
            card = card_lines
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [print_field_16(self._parse_dynamic_syntax(field)) if '%' in
                          field.strip()[0:1] else print_field_16(field) for field in fields]
                has_none = False

            if has_none:
                card = wipe_empty_fields([print_field_16(field) for field in fields])
            else:
                #card = remove_trailing_fields(fields)
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card, has_none=False)
        return card_obj, card

    def create_card_object_list(self, card_lines, card_name, has_none=True):
        """
        Creates a BDFCard object, which is really just a list that
        allows indexing past the last field

        Parameters
        ----------
        card_lines: list[str]
            the list of the card lines
            input is list of lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']
        card_name : str
            the card_name -> 'GRID'
        has_none : bool; default=True
            ???

        Returns
        -------
        card_obj : BDFCard
            the BDFCard object
        card : list[str]
            the card with empty fields removed
        """
        card_name = card_name.upper()
        self.increase_card_count(card_name)
        if card_name in ['DEQATN', 'PBRSECT', 'PBMSECT']:
            card_obj = card_lines
            card = card_lines
        else:
            fields = card_lines

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [print_field_16(self._parse_dynamic_syntax(field)) if '%' in
                          field.strip()[0:1] else print_field_16(field) for field in fields]
                has_none = False

            if has_none:
                card = wipe_empty_fields([print_field_16(field) for field in fields])
            else:
                #card = remove_trailing_fields(fields)
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card, has_none=False)
        return card_obj, card

    def create_card_object_fields(self, card_lines, card_name, has_none=True):
        """
        Creates a BDFCard object, which is really just a list that
        allows indexing past the last field

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is list of fields -> ['GRID', '1', '2', '3.0', '4.0', '5.0']
        card_name : str
            the card_name -> 'GRID'
        has_none : bool; default=True
            can there be trailing Nones in the card data
            (e.g. ['GRID', '1', '2', '3.0', '4.0', '5.0'])

        Returns
        -------
        card_obj : BDFCard
            the BDFCard object
        card : list[str]
            the card with empty fields removed
        """
        card_name = card_name.upper()
        self.increase_card_count(card_name)
        if card_name in ['DEQATN', 'PBRSECT', 'PBMSECT']:
            card_obj = card_lines
            card = card_lines
        else:
            fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [print_field_16(self._parse_dynamic_syntax(field)) if '%' in
                          field.strip()[0:1] else print_field_16(field) for field in fields]
                has_none = False

            if has_none:
                card = wipe_empty_fields([print_field_16(field) for field in fields])
            else:
                #card = remove_trailing_fields(fields)
                card = wipe_empty_fields(fields)
            card_obj = BDFCard(card, has_none=False)
        return card_obj, card

    def _make_card_parser(self):
        """creates the card parser variables that are used by add_card"""
        class Crash:
            """class for crashing on specific cards"""
            def __init__(self):
                """dummy init"""
                pass
            @classmethod
            def add_card(cls, card, comment=''):
                """the method that forces the crash"""
                raise NotImplementedError(card)

        add_methods = self #._add_methods
        self._card_parser = {
            #'=' : (Crash, None),
            '/' : (Crash, None),
            # nodes
            #'GRID' : (GRID, self.add_node),
            #'SPOINT' : (SPOINTs, self.add_spoint),
            #'EPOINT' : (EPOINTs, self.add_epoint),
            #'POINT' : (POINT, self.add_point),

            'PARAM' : (PARAM, add_methods.add_param_object),

            #'CORD2R' : (CORD2R, self._add_coord_object),
            #'CORD2C' : (CORD2C, self._add_coord_object),
            #'CORD2S' : (CORD2S, self._add_coord_object),
            #'GMCORD' : (GMCORD, self._add_coord_object),

            'PLOTEL' : (PLOTEL, add_methods.add_plotel_object),

            #'CONROD' : (CONROD, self.add_element),
            #'CROD' : (CROD, self.add_element),
            #'PROD' : (PROD, self.add_property),
            #'CTUBE' : (CTUBE, self.add_element),
            #'PTUBE' : (PTUBE, self.add_property),

            #'CBAR' : (CBAR, self.add_element),
            #'PBAR' : (PBAR, self.add_property),
            #'PBARL' : (PBARL, self.add_property),
            #'PBRSECT' : (PBRSECT, self.add_property),

            #'CBEAM' : (CBEAM, self.add_element),
            #'PBEAM' : (PBEAM, self.add_property),
            #'PBEAML' : (PBEAML, self.add_property),
            #'PBCOMP' : (PBCOMP, self.add_property),
            #'PBMSECT' : (PBMSECT, self.add_property),

            #'CBEAM3' : (CBEAM3, self.add_element),
            #'PBEAM3' : (PBEAM3, self.add_property),

            #'CBEND' : (CBEND, self.add_element),
            #'PBEND' : (PBEND, self.add_property),

            #'CTRIA3' : (CTRIA3, self.add_element),
            #'CQUAD4' : (CQUAD4, self.add_element),
            #'CQUAD' : (CQUAD, self.add_element),
            #'CQUAD8' : (CQUAD8, self.add_element),
            #'CQUADX' : (CQUADX, self.add_element),
            #'CQUADR' : (CQUADR, self.add_element),
            #'CTRIA6' : (CTRIA6, self.add_element),
            #'CTRIAR' : (CTRIAR, self.add_element),
            #'CTRIAX' : (CTRIAX, self.add_element),
            #'CTRIAX6' : (CTRIAX6, self.add_element),
            #'PCOMP' : (PCOMP, self.add_property),
            #'PCOMPG' : (PCOMPG, self.add_property),
            #'PSHELL' : (PSHELL, self.add_property),
            #'PLPLANE' : (PLPLANE, self.add_property),

            #'CPLSTN3' : (CPLSTN3, self.add_element),
            #'CPLSTN4' : (CPLSTN4, self.add_element),
            #'CPLSTN6' : (CPLSTN6, self.add_element),
            #'CPLSTN8' : (CPLSTN8, self.add_element),
            #'PPLANE' : (PPLANE, self.add_property),

            #'CSHEAR' : (CSHEAR, self.add_element),
            #'PSHEAR' : (PSHEAR, self.add_property),

            #'CTETRA' : (CTETRA, self.add_element),
            #'CPYRAM' : (CPYRAM, self.add_element),
            #'CPENTA' : (CPENTA, self.add_element),
            #'CHEXA' : (CHEXA, self.add_element),
            #'PSOLID' : (PSOLID, self.add_property),
            #'PLSOLID' : (PLSOLID, self.add_property),
            #'PCOMPS' : (PCOMPS, self.add_property),

            #'CELAS1' : (CELAS1, self.add_element),
            #'CELAS2' : (CELAS2, self.add_element),
            #'CELAS3' : (CELAS3, self.add_element),
            #'CELAS4' : (CELAS4, self.add_element),
            #'CVISC' : (CVISC, self.add_element),
            #'PELAST' : (PELAST, self.add_PELAST),

            #'CDAMP1' : (CDAMP1, self.add_damper),
            #'CDAMP2' : (CDAMP2, self.add_damper),
            #'CDAMP3' : (CDAMP3, self.add_damper),
            # CDAMP4 added later because the documentation is wrong
            #'CDAMP5' : (CDAMP5, self.add_damper),
            #'PDAMP5' : (PDAMP5, self.add_property),

            #'CFAST' : (CFAST, self.add_damper),
            #'PFAST' : (PFAST, self.add_property),

            #'CGAP' : (CGAP, self.add_element),
            #'PGAP' : (PGAP, self.add_property),

            #'CBUSH' : (CBUSH, self.add_damper),
            #'CBUSH1D' : (CBUSH1D, self.add_damper),
            #'CBUSH2D' : (CBUSH2D, self.add_damper),
            #'PBUSH' : (PBUSH, self.add_property),
            #'PBUSH1D' : (PBUSH1D, self.add_property),

            #'CRAC2D' : (CRAC2D, self.add_element),
            #'PRAC2D' : (PRAC2D, self.add_property),

            #'CRAC3D' : (CRAC3D, self.add_element),
            #'PRAC3D' : (PRAC3D, self.add_property),

            #'PDAMPT' : (PDAMPT, self.add_PDAMPT),
            #'PBUSHT' : (PBUSHT, self.add_PBUSHT),

            #'PCONEAX' : (PCONEAX, self.add_property),

            'RBAR' : (RBAR, add_methods.add_rigid_element_object),
            'RBAR1' : (RBAR1, add_methods.add_rigid_element_object),
            'RBE1' : (RBE1, add_methods.add_rigid_element_object),
            'RBE2' : (RBE2, add_methods.add_rigid_element_object),
            'RBE3' : (RBE3, add_methods.add_rigid_element_object),
            'RROD' : (RROD, add_methods.add_rigid_element_object),
            'RSPLINE' : (RSPLINE, add_methods.add_rigid_element_object),


            ## there is no MAT6 or MAT7
            #'MAT1' : (MAT1, self.add_structural_material),
            #'MAT2' : (MAT2, self.add_structural_material),
            #'MAT3' : (MAT3, self.add_structural_material),
            #'MAT8' : (MAT8, self.add_structural_material),
            #'MAT9' : (MAT9, self.add_structural_material),
            #'MAT10' : (MAT10, self.add_structural_material),
            #'MAT11' : (MAT11, self.add_structural_material),
            #'EQUIV' : (EQUIV, self.add_structural_material),

            #'MATHE' : (MATHE, self.add_hyperelastic_material),
            #'MATHP' : (MATHP, self.add_hyperelastic_material),
            #'MAT4' : (MAT4, self.add_thermal_material),
            #'MAT5' : (MAT5, self.add_thermal_material),

            #'MATS1' : (MATS1, self.add_material_dependence),
            ##'MATS3' : (MATS3, self.add_material_dependence),
            ##'MATS8' : (MATS8, self.add_material_dependence),
            #'MATT1' : (MATT1, self.add_material_dependence),
            #'MATT2' : (MATT2, self.add_material_dependence),
            ##'MATT3' : (MATT3, self.add_material_dependence),
            #'MATT4' : (MATT4, self.add_material_dependence),
            #'MATT5' : (MATT5, self.add_material_dependence),
            ##'MATT8' : (MATT8, self.add_material_dependence),
            ##'MATT9' : (MATT9, self.add_material_dependence),

            ## hasn't been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            #'CREEP' : (CREEP, self.add_creep_material),

            #'CONM1' : (CONM1, self.add_mass),
            #'CONM2' : (CONM2, self.add_mass),
            #'CMASS1' : (CMASS1, self.add_mass),
            #'CMASS2' : (CMASS2, self.add_mass),
            #'CMASS3' : (CMASS3, self.add_mass),
            ## CMASS4 - added later because documentation is wrong

            #'MPC' : (MPC, self.add_constraint_MPC),
            #'MPCADD' : (MPCADD, self.add_constraint_MPC),

            #'SPC' : (SPC, self.add_constraint_SPC),
            #'SPC1' : (SPC1, self.add_constraint_SPC1),
            #'SPCAX' : (SPCAX, self.add_constraint_SPC),
            #'SPCADD' : (SPCADD, self.add_constraint_SPC),
            #'GMSPC' : (GMSPC, self.add_constraint_SPC),

            #'SESUP' : (SESUP, self.add_sesuport), # pseudo-constraint
            #'SUPORT' : (SUPORT, self.add_suport), # pseudo-constraint
            #'SUPORT1' : (SUPORT1, self.add_suport1),  # pseudo-constraint

            #'FORCE' : (FORCE, self.add_load),
            #'FORCE1' : (FORCE1, self.add_load),
            #'FORCE2' : (FORCE2, self.add_load),
            #'MOMENT' : (MOMENT, self.add_load),
            #'MOMENT1' : (MOMENT1, self.add_load),
            #'MOMENT2' : (MOMENT2, self.add_load),

            #'LSEQ' : (LSEQ, self.add_LSEQ),
            #'LOAD' : (LOAD, self.add_load),
            #'LOADCYN' : (LOADCYN, self.add_load),
            #'GRAV' : (GRAV, self.add_load),
            #'ACCEL' : (ACCEL, self.add_load),
            #'ACCEL1' : (ACCEL1, self.add_load),
            #'PLOAD' : (PLOAD, self.add_load),
            #'PLOAD1' : (PLOAD1, self.add_load),
            #'PLOAD2' : (PLOAD2, self.add_load),
            #'PLOAD4' : (PLOAD4, self.add_load),
            #'PLOADX1' : (PLOADX1, self.add_load),
            #'RFORCE' : (RFORCE, self.add_load),
            #'RFORCE1' : (RFORCE1, self.add_load),
            #'SLOAD' : (SLOAD, self.add_load),
            #'RANDPS' : (RANDPS, self.add_load),
            #'GMLOAD' : (GMLOAD, self.add_load),
            #'SPCD' : (SPCD, self.add_load),  # enforced displacement
            #'QVOL' : (QVOL, self.add_load),  # thermal

            #'DLOAD' : (DLOAD, self.add_dload),
            #'ACSRCE' : (ACSRCE, self._add_dload_entry),
            #'TLOAD1' : (TLOAD1, self._add_dload_entry),
            #'TLOAD2' : (TLOAD2, self._add_dload_entry),
            #'RLOAD1' : (RLOAD1, self._add_dload_entry),
            #'RLOAD2' : (RLOAD2, self._add_dload_entry),

            #'FREQ' : (FREQ, self.add_FREQ),
            #'FREQ1' : (FREQ1, self.add_FREQ),
            #'FREQ2' : (FREQ2, self.add_FREQ),
            #'FREQ4' : (FREQ4, self.add_FREQ),

            'DOPTPRM' : (DOPTPRM, add_methods.add_doptprm_object),
            'DESVAR' : (DESVAR, add_methods.add_desvar_object),
            # BCTSET

            #'TEMP' : (TEMP, self.add_thermal_load),
            #'QBDY1' : (QBDY1, self.add_thermal_load),
            #'QBDY2' : (QBDY2, self.add_thermal_load),
            #'QBDY3' : (QBDY3, self.add_thermal_load),
            #'QHBDY' : (QHBDY, self.add_thermal_load),
            #'PHBDY' : (PHBDY, self.add_PHBDY),

            #'CHBDYE' : (CHBDYE, self.add_thermal_element),
            #'CHBDYG' : (CHBDYG, self.add_thermal_element),
            #'CHBDYP' : (CHBDYP, self.add_thermal_element),
            #'PCONV' : (PCONV, self.add_convection_property),
            #'PCONVM' : (PCONVM, self.add_convection_property),

            # aero
            'AECOMP' : (AECOMP, add_methods.add_aecomp_object),
            'AEFACT' : (AEFACT, add_methods.add_aefact_object),
            'AELINK' : (AELINK, add_methods.add_aelink_object),
            'AELIST' : (AELIST, add_methods.add_aelist_object),
            'AEPARM' : (AEPARM, add_methods.add_aeparm_object),
            'AESTAT' : (AESTAT, add_methods.add_aestat_object),
            'AESURF' : (AESURF, add_methods.add_aesurf_object),
            'AESURFS' : (AESURFS, add_methods.add_aesurfs_object),

            'CAERO1' : (CAERO1, add_methods.add_caero_object),
            'CAERO2' : (CAERO2, add_methods.add_caero_object),
            'CAERO3' : (CAERO3, add_methods.add_caero_object),
            'CAERO4' : (CAERO4, add_methods.add_caero_object),
            'CAERO5' : (CAERO5, add_methods.add_caero_object),

            'PAERO1' : (PAERO1, add_methods.add_paero_object),
            'PAERO2' : (PAERO2, add_methods.add_paero_object),
            'PAERO3' : (PAERO3, add_methods.add_paero_object),
            'PAERO4' : (PAERO4, add_methods.add_paero_object),
            'PAERO5' : (PAERO5, add_methods.add_paero_object),

            'SPLINE1' : (SPLINE1, add_methods.add_spline_object),
            'SPLINE2' : (SPLINE2, add_methods.add_spline_object),
            'SPLINE3' : (SPLINE3, add_methods.add_spline_object),
            'SPLINE4' : (SPLINE4, add_methods.add_spline_object),
            'SPLINE5' : (SPLINE5, add_methods.add_spline_object),

            # SOL 144
            'AEROS' : (AEROS, add_methods.add_aeros_object),
            'TRIM' : (TRIM, add_methods.add_trim_object),
            'DIVERG' : (DIVERG, add_methods.add_diverg_object),

            # SOL 145
            'AERO' : (AERO, add_methods.add_aero_object),
            'FLUTTER' : (FLUTTER, add_methods.add_flutter_object),
            'FLFACT' : (FLFACT, add_methods.add_flfact_object),
            'MKAERO1' : (MKAERO1, add_methods.add_mkaero_object),
            'MKAERO2' : (MKAERO2, add_methods.add_mkaero_object),

            'GUST' : (GUST, add_methods.add_gust_object),
            'CSSCHD' : (CSSCHD, add_methods.add_csschd_object),
            'MONPNT1' : (MONPNT1, add_methods.add_monpnt_object),
            'MONPNT2' : (MONPNT2, add_methods.add_monpnt_object),
            'MONPNT3' : (MONPNT3, add_methods.add_monpnt_object),

            'NLPARM' : (NLPARM, add_methods.add_nlparm_object),
            'NLPCI' : (NLPCI, add_methods.add_nlpci_object),
            'TSTEP' : (TSTEP, add_methods.add_tstep_object),
            'TSTEPNL' : (TSTEPNL, add_methods.add_tstepnl_object),

            #'TF' : (TF, self.add_TF),
            #'DELAY' : (DELAY, self.add_DELAY),

            'DCONADD' : (DCONADD, add_methods.add_dconstr_object),
            'DCONSTR' : (DCONSTR, add_methods.add_dconstr_object),
            'DDVAL' : (DDVAL, add_methods.add_ddval_object),
            'DLINK' : (DLINK, add_methods.add_dlink_object),

            #'DTABLE' : (DTABLE, self.add_dtable),
            'DRESP1' : (DRESP1, add_methods.add_dresp_object),
            'DRESP2' : (DRESP2, add_methods.add_dresp_object), # deqatn
            'DRESP3' : (DRESP3, add_methods.add_dresp_object),
            'DVCREL1' : (DVCREL1, add_methods.add_dvcrel_object), # dvcrels
            'DVCREL2' : (DVCREL2, add_methods.add_dvcrel_object),
            'DVPREL1' : (DVPREL1, add_methods.add_dvprel_object), # dvprels
            'DVPREL2' : (DVPREL2, add_methods.add_dvprel_object),
            'DVMREL1' : (DVMREL1, add_methods.add_dvmrel_object), # ddvmrels
            'DVMREL2' : (DVMREL2, add_methods.add_dvmrel_object),
            'DVGRID' : (DVGRID, add_methods.add_dvgrid_object), # dvgrids

            #'TABLED1' : (TABLED1, self.add_table),
            #'TABLED2' : (TABLED2, self.add_table),
            #'TABLED3' : (TABLED3, self.add_table),
            #'TABLED4' : (TABLED4, self.add_table),
            #'TABLEM1' : (TABLEM1, self.add_table),
            #'TABLEM2' : (TABLEM2, self.add_table),
            #'TABLEM3' : (TABLEM3, self.add_table),
            #'TABLEM4' : (TABLEM4, self.add_table),

            #'TABLES1' : (TABLES1, self.add_table),
            #'TABLEST' : (TABLEST, self.add_table),

            #'TABDMP1' : (TABDMP1, self.add_table_sdamping),
            #'TABRND1' : (TABRND1, self.add_random_table),
            #'TABRNDG' : (TABRNDG, self.add_random_table),

            'EIGB' : (EIGB, add_methods.add_method_object),
            'EIGR' : (EIGR, add_methods.add_method_object),
            'EIGRL' : (EIGRL, add_methods.add_method_object),
            'EIGC' : (EIGC, add_methods.add_cmethod_object),
            'EIGP' : (EIGP, add_methods.add_cmethod_object),

            'BCRPARA' : (BCRPARA, add_methods.add_bcrpara_object),
            'BCTADD' : (BCTADD, add_methods.add_bctadd_object),
            'BCTPARA' : (BCTPARA, add_methods.add_bctpara_object),
            'BSURF' : (BSURF, add_methods.add_bsurf_object),
            'BSURFS' : (BSURFS, add_methods.add_bsurfs_object),

            'ASET' : (ASET, add_methods.add_aset_object),
            'ASET1' : (ASET1, add_methods.add_aset_object),

            'BSET' : (BSET, add_methods.add_bset_object),
            'BSET1' : (BSET1, add_methods.add_bset_object),

            'CSET' : (CSET, add_methods.add_cset_object),
            'CSET1' : (CSET1, add_methods.add_cset_object),

            'QSET' : (QSET, add_methods.add_qset_object),
            'QSET1' : (QSET1, add_methods.add_qset_object),

            'USET' : (USET, add_methods.add_uset_object),
            'USET1' : (USET1, add_methods.add_uset_object),

            'SET1' : (SET1, add_methods.add_set_object),
            'SET3' : (SET3, add_methods.add_set_object),

            'SESET' : (SESET, add_methods.add_seset_object),

            'SEBSET' : (SEBSET, add_methods.add_sebset_object),
            'SEBSET1' : (SEBSET1, add_methods.add_sebset_object),

            'SECSET' : (SECSET, add_methods.add_secset_object),
            'SECSET1' : (SECSET1, add_methods.add_secset_object),

            'SEQSET' : (SEQSET, add_methods.add_seqset_object),
            'SEQSET1' : (SEQSET1, add_methods.add_seqset_object),

            #'SESUP' : (SESUP, self.add_SESUP),  # pseudo-constraint

            #'SEUSET' : (SEUSET, self.add_SEUSET),
            #'SEUSET1' : (SEUSET1, self.add_SEUSET),

            # BCTSET
        }
        self._card_parser_prepare = {
            #'CORD2R' : (CORD2R, self._add_coord_object), # not vectorized
            #'CORD2C' : (CORD2C, self._add_coord_object),
            #'CORD2S' : (CORD2S, self._add_coord_object),
            'CORD2R' : self._prepare_cord2, # vectorized
            'CORD2C' : self._prepare_cord2,
            'CORD2S' : self._prepare_cord2,

            #'CORD1R' : self._prepare_cord1r,
            #'CORD1C' : self._prepare_cord1c,
            #'CORD1S' : self._prepare_cord1s,
            ##'CORD3G' : self._prepare_CORD3G,

            #'DAREA' : self._prepare_darea,
            #'DPHASE' : self._prepare_dphase,
            #'PMASS' : self._prepare_pmass,
            #'CMASS4' : self._prepare_cmass4,
            #'CDAMP4' : self._prepare_cdamp4,

            'DMIG' : self._prepare_dmig,
            'DMI' : self._prepare_dmi,
            'DMIJ' : self._prepare_dmij,
            'DMIK' : self._prepare_dmik,
            'DMIJI' : self._prepare_dmiji,

            'DEQATN' : self._prepare_dequatn,

            #'PVISC' : self._prepare_pvisc,
            #'PELAS' : self._prepare_pelas,
            #'PDAMP' : self._prepare_pdamp,

            #'TEMPD' : self._prepare_tempd,
            #'CONVM' : self._prepare_convm,
            #'CONV' : self._prepare_conv,
            #'RADM' : self._prepare_radm,
            #'RADBC' : self._prepare_radbc,
            ## GRDSET-will be last card to update from _card_parser_prepare
            #'GRDSET' : self._prepare_grdset,

            #'BCTSET' : self._prepare_bctset,
        }

    def reject_card_obj2(self, card_name, card_obj):
        """rejects a card object"""
        self.reject_cards.append(card_obj)

    def reject_card_lines(self, card_name: str, card_lines: list[str],
                          show_log: bool=True, comment: str='') -> None:
        """rejects a card"""
        if card_name.isdigit():
            # TODO: this should technically work (I think), but it's a problem
            #       for the code
            #
            # prevents:
            # spc1,100,456,10013832,10013833,10013830,10013831,10013836,10013837,
            # 10013834,10013835,10013838,10013839,10014508,10008937,10008936,10008935,
            msg = 'card_name=%r was misparsed...\ncard_lines=%s' % (
                card_name, card_lines)
            raise RuntimeError(msg)
        if card_name not in self.card_count:
            if ' ' in card_name:
                msg = (
                    'No spaces allowed in card name %r.  '
                    'Should this be a comment?\n%s%s' % (
                        card_name, comment, card_lines))
                raise RuntimeError(msg)
            if card_name in ['SUBCASE ', 'CEND']:
                raise RuntimeError('No executive/case control deck was defined.')
            self.log.info(f'    rejecting card_name = {card_name}')
            if len(card_name) <= 3:
                self.log.error(f'    card_lines = {card_lines}')
        self.increase_card_count(card_name)
        self.rejects.append([comment] + card_lines)

    def _prepare_bctset(self, card, card_obj, comment=''):
        """adds a GRDSET"""
        card = BCTSET.add_card(card_obj, comment=comment, sol=self.sol)
        self._add_bctset_object(card)

    def _prepare_grdset(self, card, card_obj, comment=''):
        """adds a GRDSET"""
        self.grdset = GRDSET.add_card(card_obj, comment=comment)

    #def _prepare_cdamp4(self, card, card_obj, comment=''):
        #"""adds a CDAMP4"""
        #self.add_damper(CDAMP4.add_card(card_obj, comment=comment))
        #if card_obj.field(5):
            #self.add_damper(CDAMP4.add_card(card_obj, 1, comment=''))
        #return card_obj

    def _prepare_convm(self, card, card_obj, comment=''):
        """adds a CONVM"""
        boundary_condition = CONVM.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)

    def _prepare_conv(self, card, card_obj, comment=''):
        """adds a CONV"""
        boundary_condition = CONV.add_card(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.eid)

    def _prepare_radm(self, card, card_obj, comment=''):
        """adds a RADM"""
        boundary_condition = RADM.add_card(card, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.radmid)

    def _prepare_radbc(self, card, card_obj, comment=''):
        """adds a RADBC"""
        boundary_condition = RADBC(card_obj, comment=comment)
        self._add_thermal_bc_object(boundary_condition, boundary_condition.nodamb)

    def _prepare_tempd(self, card, card_obj, comment=''):
        """adds a TEMPD"""
        self.add_tempd(TEMPD.add_card(card_obj, 0, comment=comment))
        if card_obj.field(3):
            self.add_tempd(TEMPD.add_card(card_obj, 1, comment=''))
            if card_obj.field(5):
                self.add_tempd(TEMPD.add_card(card_obj, 2, comment=''))
                if card_obj.field(7):
                    self.add_tempd(TEMPD.add_card(card_obj, 3, comment=''))

    def _add_doptprm(self, doptprm, comment=''):
        """adds a DOPTPRM"""
        self.doptprm = doptprm

    def _prepare_dequatn(self, card, card_obj, comment=''):
        """adds a DEQATN"""
        if hasattr(self, 'test_deqatn') or 1:
            self.add_deqatn(DEQATN.add_card(card_obj, comment=comment))
        else:
            if comment:
                self.rejects.append([comment])
            self.rejects.append(card)

    def _prepare_dmig(self, card, card_obj, comment=''):
        """adds a DMIG"""
        name = string(card_obj, 1, 'name')
        field2 = integer_or_string(card_obj, 2, 'flag')
        #print('name=%r field2=%r' % (name, field2))

        if name == 'UACCEL':  # special DMIG card
            if field2 == 0:
                card = DMIG_UACCEL.add_card(card_obj, comment=comment)
                self.add_dmig(card)
            else:
                self._dmig_temp[name].append((card_obj, comment))
        else:
            field2 = integer_or_string(card_obj, 2, 'flag')
            if field2 == 0:
                card = DMIG(card_obj, comment=comment)
                self.add_dmig(card)
            else:
                self._dmig_temp[name].append((card_obj, comment))

    def _prepare_dmix(self, class_obj, add_method, card_obj, comment=''):
        """adds a DMIx"""
        #elif card_name in ['DMI', 'DMIJ', 'DMIJI', 'DMIK']:
        field2 = integer(card_obj, 2, 'flag')
        if field2 == 0:
            add_method(class_obj(card_obj, comment=comment))
        else:
            name = string(card_obj, 1, 'name')
            self._dmig_temp[name].append((card_obj, comment))

    def _prepare_dmi(self, card, card_obj, comment=''):
        """adds a DMI"""
        self._prepare_dmix(DMI, self._add_dmi_object, card_obj, comment=comment)

    def _prepare_dmij(self, card, card_obj, comment=''):
        """adds a DMIJ"""
        self._prepare_dmix(DMIJ, self._add_dmij_object, card_obj, comment=comment)

    def _prepare_dmik(self, card, card_obj, comment=''):
        """adds a DMIK"""
        self._prepare_dmix(DMIK, self._add_dmik_object, card_obj, comment=comment)

    def _prepare_dmiji(self, card, card_obj, comment=''):
        """adds a DMIJI"""
        self._prepare_dmix(DMIJI, self._add_dmiji_object, card_obj, comment=comment)

    #def _prepare_cmass4(self, card, card_obj, comment=''):
        #"""adds a CMASS4"""
        #class_instance = CMASS4.add_card(card_obj, icard=0, comment=comment)
        #self.add_mass(class_instance)
        #if card_obj.field(5):
            #class_instance = CMASS4.add_card(card_obj, icard=1, comment=comment)
            #self.add_mass(class_instance)

    #def _prepare_pelas(self, card, card_obj, comment=''):
        #"""adds a PELAS"""
        #class_instance = PELAS.add_card(card_obj, icard=0, comment=comment)
        #self.add_property(class_instance)
        #if card_obj.field(5):
            #class_instance = PELAS.add_card(card_obj, icard=1, comment=comment)
            #self.add_property(class_instance)

    #def _prepare_pvisc(self, card, card_obj, comment=''):
        #"""adds a PVISC"""
        #class_instance = PVISC.add_card(card_obj, icard=0, comment=comment)
        #self.add_property(class_instance)
        #if card_obj.field(5):
            #class_instance = PVISC.add_card(card_obj, icard=1, comment=comment)
            #self.add_property(class_instance)

    #def _prepare_pdamp(self, card, card_obj, comment=''):
        #"""adds a PDAMP"""
        #class_instance = PDAMP.add_card(card_obj, icard=0, comment=comment)
        #self.add_property(class_instance)
        #if card_obj.field(3):
            #class_instance = PDAMP.add_card(card_obj, icard=1, comment=comment)
            #self.add_property(class_instance)
        #if card_obj.field(5):
            #class_instance = PDAMP.add_card(card_obj, icard=2, comment=comment)
            #self.add_property(class_instance)
        #if card_obj.field(7):
            #class_instance = PDAMP.add_card(card_obj, icard=3, comment=comment)
            #self.add_property(class_instance)

    #def _prepare_pmass(self, card, card_obj, comment=''):
        #"""adds a PMASS"""
        #card_instance = PMASS(card_obj, icard=0, comment=comment)
        #self.add_property_mass(card_instance)
        #for (i, j) in enumerate([3, 5, 7]):
            #if card_obj.field(j):
                #card_instance = PMASS(card_obj, icard=i+1, comment=comment)
                #self.add_property_mass(card_instance)

    #def _prepare_dphase(self, card, card_obj, comment=''):
        #"""adds a DPHASE"""
        #class_instance = DPHASE.add_card(card_obj, comment=comment)
        #self.add_dphase(class_instance)
        #if card_obj.field(5):
            #print('card_obj = ', card_obj)
            #class_instance = DPHASE(card_obj, icard=1, comment=comment)
            #self.add_DPHASE(class_instance)

    def _prepare_cord1r(self, card, card_obj, comment=''):
        """adds a CORD1R"""
        class_instance = CORD1R.add_card(card_obj, comment=comment)
        self._add_methods.add_coord_object(class_instance)
        if card_obj.field(5):
            class_instance = CORD1R.add_card(card_obj, icard=1, comment=comment)
            self._add_methods.add_coord_object(class_instance)

    def _prepare_cord1c(self, card, card_obj, comment=''):
        """adds a CORD1C"""
        class_instance = CORD1C.add_card(card_obj, comment=comment)
        self._add_methods.add_coord_object(class_instance)
        if card_obj.field(5):
            class_instance = CORD1C.add_card(card_obj, icard=1, comment=comment)
            self._add_methods.add_coord_object(class_instance)

    def _prepare_cord1s(self, card, card_obj, comment=''):
        """adds a CORD1S"""
        class_instance = CORD1S.add_card(card_obj, comment=comment)
        self._add_methods.add_coord_object(class_instance)
        if card_obj.field(5):
            class_instance = CORD1S.add_card(card_obj, icard=1, comment=comment)
            self._add_methods.add_coord_object(class_instance)

    def _prepare_cord2(self, card, card_obj, comment=''):
        """adds a CORD2x"""
        self.coords.add_cord2x(card, card_obj, comment)

    def add_card(self, card_lines, card_name, comment='', is_list=True, has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        is_list : bool, optional
            False : input is a list of card fields -> ['GRID', 1, None, 3.0, 4.0, 5.0]
            True :  input is list of card_lines -> ['GRID, 1,, 3.0, 4.0, 5.0']
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card

        .. code-block:: python

           >>> model = BDF()

           # is_list is a somewhat misleading name; is it a list of card_lines
           # where a card_line is an unparsed string
           >>> card_lines = ['GRID,1,2']
           >>> comment = 'this is a comment'
           >>> model.add_card(card_lines, 'GRID', comment, is_list=True)

           # here is_list=False because it's been parsed
           >>> card = ['GRID', 1, 2,]
           >>> model.add_card(card_lines, 'GRID', comment, is_list=False)

           # here is_list=False because it's been parsed
           # Note the None at the end of the 1st line, which is there
           #      because the CONM2 card has a blank field.
           #      It must be there.
           # We also set i32 on the 2nd line, so it will default to 0.0
           >>> card = [
                   'CONM2', eid, nid, cid, mass, x1, x2, x3, None,
                            i11, i21, i22, i31, None, i33,
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=False)

           # here's an alternate approach for the CONM2
           # we use Nastran's CSV format
           # There are many blank fields, but it's parsed exactly like a
           # standard CONM2.
           >>> card = [
                   'CONM2,1,2,3,10.0',
                   ',1.0,,5.0'
               ]
           >>> model.add_card(card_lines, 'CONM2', comment, is_list=True)

        Notes
        -----
        This is a very useful method for interfacing with the code.
        The card_object is not a card-type object...so not a GRID
        card or CQUAD4 object.  It's a BDFCard Object.  However,
        you know the type (assuming a GRID), so just call the
        *mesh.Node(nid)* to get the Node object that was just
        created.
        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name,
                                                 is_list=is_list, has_none=has_none)
        self._add_card_helper(card_obj, card_name, card_name, comment)
        return card_obj

    def add_card_fields(self, card_lines, card_name, comment='', has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is a list of card fields -> ['GRID', 1, 2, 3.0, 4.0, 5.0]
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID', 1, 2, 3.0, 4.0, 5.0, None])

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card
        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name,
                                                 is_list=True, has_none=has_none)
        self._add_card_helper(card_obj, card, card_name, comment)

    @property
    def nodes(self):
        ngrids = len(self.grid)
        assert ngrids > 0, ngrids
        nspoints = 0
        nepoints = 0
        if self.spoint.n:
            spoints = self.spoint.points
            nspoints = len(spoints)
        if self.epoint.n:
            epoints = self.epoint.points
            nepoints = len(epoints)
            raise NotImplementedError('EPOINT')

        assert ngrids + nspoints + nepoints > 0, 'ngrids=%s nspoints=%s nepoints=%s' % (ngrids, nspoints, nepoints)
        nodes = np.zeros(ngrids + nspoints + nepoints, dtype='int32')
        nodes[:ngrids] = self.grid.node_id
        if nspoints:
            nodes[ngrids:ngrids+nspoints] = self.spoint.points
        if nepoints:
            nodes[ngrids+nspoints:] = self.epoint.points
        return nodes

    def get_xyz_in_coord(self, cid=0, fdtype='float64', sort_ids=True):
        """
        Gets the xyz points (including SPOINTS) in the desired coordinate frame

        Parameters
        ----------
        cid : int; default=0
            the desired coordinate system
        fdtype : str; default='float64'
            the data type of the xyz coordinates
        sort_ids : bool; default=True
            sort the ids

        Returns
        -------
        xyz : (n, 3) ndarray
            the xyz points in the cid coordinate frame

        .. warning:: doesn't support EPOINTs
        """
        ngrids = len(self.grid)
        nspoints = 0
        nepoints = 0
        spoints = None
        if self.spoint.n:
            spoints = self.point.points
            nspoints = len(spoints)
        if self.epoint.n:
            epoints = self.point.points
            nepoints = len(epoints)
            raise NotImplementedError('EPOINT')

        assert ngrids + nspoints + nepoints > 0, 'ngrids=%s nspoints=%s nepoints=%s' % (ngrids, nspoints, nepoints)
        xyz_cid0 = np.zeros((ngrids + nspoints + nepoints, 3), dtype=fdtype)
        if cid == 0:
            xyz_cid0 = self.grid.get_position_by_node_index()
            assert nspoints == 0, nspoints
        else:
            assert cid == 0, cid
            assert nspoints == 0, nspoints
        if sort_ids:
            isort = np.argsort(all_nodes)
            xyz_cid0 = xyz_cid0[isort, :]

        return xyz_cid0

    def _add_card_helper(self, card_obj, card, card_name, comment=''):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_object : BDFCard()
            the card object representation of card
        card : list[str]
            the fields of the card object; used for rejection and special cards
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        """
        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo:
            try:
                print(print_card_8(card_obj).rstrip())
            except Exception:
                print(print_card_16(card_obj).rstrip())

        if card_name in self._card_parser:
            card_class, add_card_function = self._card_parser[card_name]
            try:
                class_instance = card_class.add_card(card_obj, comment=comment)
                add_card_function(class_instance)
            except TypeError:
                #msg = 'problem adding %s' % card_obj
                raise
                #raise TypeError(msg)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                #self._iparse_errors += 1
                #self.log.error(card_obj)
                #var = traceback.format_exception_only(type(exception), exception)
                #self._stored_parse_errors.append((card, var))
                #if self._iparse_errors > self._nparse_errors:
                    #self.pop_parse_errors()
                #raise
            #except AssertionError as exception:
                #self.log.error(card_obj)

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            try:
                add_card_function(card, card_obj, comment)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                #self._iparse_errors += 1
                #self.log.error(card_obj)
                #var = traceback.format_exception_only(type(exception), exception)
                #self._stored_parse_errors.append((card, var))
                #if self._iparse_errors > self._nparse_errors:
                    #self.pop_parse_errors()
            #except AssertionError as exception:
                #self.log.error(card_obj)
                #raise
        else:
            #raise RuntimeError(card_obj)
            self.reject_cards.append(card_obj)

    def get_bdf_stats(self, return_type='string'):
        """
        Print statistics for the BDF

        Parameters
        ----------
        return_type : str (default='string')
            the output type ('list', 'string')
                'list' : list of strings
                'string' : single, joined string

        Returns
        -------
        return_data : str, optional
            the output data

            Notes
        -----
        If a card is not supported and not added to the proper
        lists, this method will fail.
        """
        return ''
        card_stats = [
            'params', 'nodes', 'points', 'elements', 'rigid_elements',
            'properties', 'materials', 'creep_materials',
            'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1', 'MATS3', 'MATT8',
            'coords', 'mpcs', 'mpcadds',

            # dynamic cards
            'dareas', 'dphases', 'nlparms', 'nlpcis', 'tsteps', 'tstepnls',

            # direct matrix input - DMIG - dict
            'dmi', 'dmig', 'dmij', 'dmiji', 'dmik',
            'dequations',

            # frequencies - dict
            'frequencies',

            # optimization - dict
            'dconadds', 'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps',
            'dvcrels', 'dvmrels', 'dvprels', 'dvgrids',

            # SESETx - dict
            'suport1',
            'se_sets',
            'se_usets',

            # tables
            'tables', 'tables_d', 'tables_m', 'random_tables',

            # methods
            'methods', 'cMethods',

            # aero
            'caeros', 'paeros', 'aecomps', 'aefacts', 'aelinks',
            'aelists', 'aeparams', 'aesurfs', 'aestats', 'gusts', 'flfacts',
            'flutters', 'splines', 'trims',

            # thermal
            'bcs', 'thermal_materials', 'phbdys',
            'convection_properties', ]

        # These are ignored because they're lists
        ignored_types = set([
            'spoints', 'spointi',  # singleton
            'grdset',  # singleton

            'spcs',

            'suport', 'se_suport', # suport, suport1 - list
            'doptprm',  # singleton

            # SETx - list
            'sets', 'asets', 'bsets', 'csets', 'qsets',
            'se_bsets', 'se_csets', 'se_qsets',
        ])

        ## TODO: why are some of these ignored?
        ignored_types2 = set([
            'case_control_deck', 'caseControlDeck',
            'spcObject2', 'mpcObject2',

            # done
            'sol', 'loads', 'mkaeros',
            'rejects', 'reject_cards',

            # not cards
            'debug', 'executive_control_lines',
            'case_control_lines', 'cards_to_read', 'card_count',
            'is_structured', 'uniqueBulkDataCards',
            'model_type', 'include_dir',
            'sol_method', 'log',
            'sol_iline',
            'reject_count', '_relpath',
            #'foundEndData',
            'special_cards',])

        unsupported_types = ignored_types.union(ignored_types2)
        all_params = object_attributes(self, keys_to_skip=unsupported_types)

        msg = [
            '---BDF Statistics---',
            'SOL %s\n' % self.sol,
        ]

        # loads
        for (lid, loads) in sorted(self.loads.items()):
            msg.append('bdf.loads[%s]' % lid)
            groups_dict = {}
            for loadi in loads:
                groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        # dloads
        for (lid, loads) in sorted(self.dloads.items()):
            msg.append('bdf.dloads[%s]' % lid)
            groups_dict = {}
            for loadi in loads:
                groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        for (lid, loads) in sorted(self.dload_entries.items()):
            msg.append('bdf.dload_entries[%s]' % lid)
            groups_dict = {}
            for loadi in loads:
                groups_dict[loadi.type] = groups_dict.get(loadi.type, 0) + 1
            for name, count_name in sorted(groups_dict.items()):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        # aero
        if self.aero:
            msg.append('bdf:aero')
            msg.append('  %-8s %s' % ('AERO:', 1))

        # aeros
        if self.aeros:
            msg.append('bdf:aeros')
            msg.append('  %-8s %s' % ('AEROS:', 1))

        #mkaeros
        if self.mkaeros:
            msg.append('bdf:mkaeros')
            msg.append('  %-8s %s' % ('MKAERO:', len(self.mkaeros)))

        for card_group_name in card_stats:
            card_group = getattr(self, card_group_name)
            groups = set()

            if not isinstance(card_group, dict):
                msg = '%s is a %s; not dictionary' % (card_group_name, type(card_group))
                raise RuntimeError(msg)
            for card in card_group.values():
                if isinstance(card, list):
                    for card2 in card:
                        groups.add(card2.type)
                else:
                    groups.add(card.type)

            group_msg = []
            for card_name in sorted(groups):
                try:
                    ncards = self.card_count[card_name]
                    group_msg.append('  %-8s : %s' % (card_name, ncards))
                except KeyError:
                    group_msg.append('  %-8s : ???' % card_name)
                    #assert card_name == 'CORD2R', self.card_count
            if group_msg:
                msg.append('bdf.%s' % card_group_name)
                msg.append('\n'.join(group_msg))
                msg.append('')

        # rejects
        if self.rejects:
            msg.append('Rejected Cards')
            for name, counter in sorted(self.card_count.items()):
                if name not in self.cards_to_read:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        if return_type == 'string':
            return '\n'.join(msg)
        else:
            return msg

    def get_displacement_index_xyz_cp_cd(self, fdtype='float64', idtype='int32'):
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Parameters
        ----------
        fdtype : str; default='float64'
            the type of xyz_cp
        idtype : str; default='int32'
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

        Examples
        --------
        # assume GRID 1 has a CD=10
        # assume GRID 2 has a CD=10
        # assume GRID 5 has a CD=50
        >>> model.point_ids
        [1, 2, 5]
        >>> i_transform = model.get_displacement_index_xyz_cp_cd()
        >>> i_transform[10]
        [0, 1]

        >>> i_transform[50]
        [2]
        """
        nids_cd_transform = defaultdict(list)
        nids_cp_transform = defaultdict(list)
        i_transform = {}

        nnodes = len(self.nodes)
        nspoints = 0
        nepoints = 0
        spoints = None
        epoints = None
        if 0 and self.new_spoints:
            if self.new_spoints:
                if self.spoints:
                    spoints = list(self.spoints)
                    nspoints = len(spoints)
                    all_nodes += spoints
                if self.epoints:
                    epoints = list(self.epoints)
                    nepoints = len(epoints)
                    all_nodes += epoints
        else:
            if self.spoints:
                spoints = self.spoints.points
                nspoints = len(spoints)
            if self.epoints is not None:
                epoints = self.epoints.points
                nepoints = len(epoints)
                #raise NotImplementedError('EPOINTs')

        if nnodes + nspoints + nepoints == 0:
            msg = f'nnodes={nnodes:d} nspoints={nspoints:d} nepoints={nepoints:d}'
            raise ValueError(msg)

        #xyz_cid0 = np.zeros((nnodes + nspoints, 3), dtype=dtype)
        xyz_cp = np.zeros((nnodes + nspoints, 3), dtype=fdtype)
        nid_cp_cd = np.zeros((nnodes + nspoints, 3), dtype=idtype)
        i = 0
        for nid, node in sorted(self.nodes.items()):
            cd = node.Cd()
            cp = node.Cp()
            nids_cd_transform[cp].append(nid)
            nids_cd_transform[cd].append(nid)
            nid_cp_cd[i, :] = [nid, cp, cd]
            xyz_cp[i, :] = node.xyz
            i += 1
        if nspoints:
            for nid in sorted(self.spoints.points):
                nid_cp_cd[i] = nid
                i += 1
        if nepoints:
            for nid in sorted(self.epoints.points):
                nid_cp_cd[i] = nid
                i += 1

        if sort_ids:
            nids = nid_cp_cd[:, 0]
            isort = nids.argsort()
            nid_cp_cd = nid_cp_cd[isort, :]
            xyz_cp = xyz_cp[isort, :]

        icp_transform = {}
        icd_transform = {}
        nids_all = np.array(sorted(self.point_ids))
        for cd, nids in sorted(nids_cd_transform.items()):
            if cd in -1:
                continue
            nids = np.array(nids)
            icd_transform[cd] = np.where(np.isin(nids_all, nids))[0]
            if cd in nids_cp_transform:
                icp_transform[cd] = icd_transform[cd]

        for cp, nids in sorted(nids_cd_transform.items()):
            if cp in -1:
                continue
            # if cp in icd_transform:
            #     continue
            nids = np.array(nids)
            icd_transform[cp] = np.where(np.isin(nids_all, nids))[0]

        return icd_transform, icp_transform, xyz_cp, nid_cp_cd

    def transform_xyzcp_to_xyz_cid(self, xyz_cp, nids, icp_transform, in_place=False, cid=0):
        """
        Working on faster method for calculating node locations
        Not validated...

        Parameters
        ----------
        xyz_cp : (n, 3) float ndarray
            points in the CP coordinate system
        icp_transform : dict{int cp : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their input (`CP`) in that
            coordinate system.
        cid : int; default=0
            the coordinate system to get xyz in

        Returns
        -------
        xyz_cid : (n, 3) float ndarray
            points in the CID coordinate system
        """
        coord2 = self.coords[cid]
        beta2 = coord2.beta()

        assert in_place is False, 'in_place=%s' % in_place
        if in_place:
            xyz_cid0 = xyz_cp
        else:
            xyz_cid0 = np.copy(xyz_cp)

        # transform the grids to the global coordinate system
        for cp, inode in icp_transform.items():
            if cp == 0:
                continue
            coord = self.coords[cp]
            beta = coord.beta()
            is_beta = np.abs(np.diagonal(beta)).min() == 1.
            is_origin = np.abs(coord.origin).max() == 0.
            xyzi = coord.coord_to_xyz_array(xyz_cp[inode, :])
            if is_beta and is_origin:
                xyz_cid0[inode, :] = xyzi @ beta + coord.origin
            elif is_beta:
                xyz_cid0[inode, :] = xyzi @ beta
            else:
                xyz_cid0[inode, :] = xyzi + coord.origin

        if cid == 0:
            return xyz_cid0

        is_beta = np.abs(np.diagonal(beta2)).min() == 1.
        is_origin = np.abs(coord2.origin).max() == 0.
        if is_beta and is_origin:
            xyzi = (xyz_cid0 - coord2.origin) @ beta2.T
            xyz_cid = coord2.xyz_to_coord_array(xyzi)
        elif is_beta:
            xyzi = xyz_cid0 @ beta2.T
            xyz_cid = coord2.xyz_to_coord_array(xyzi)
        else:
            xyzi = xyz_cid0 - coord2.origin
            xyz_cid = coord2.xyz_to_coord_array(xyzi)

        return xyz_cid

    def get_displacement_index(self):
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Returns
        -------
        icd_transform : dict{int cid : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.

        Examples
        --------
        # assume GRID 1 has a CD=10
        # assume GRID 2 has a CD=10
        # assume GRID 5 has a CD=50
        >>> model.point_ids
        [1, 2, 5]
        >>> icd_transform = model.get_displacement_index()
        >>> icd_transform[10]
        [0, 1]

        >>> icd_transform[50]
        [2]
        """
        nids_transform = defaultdict(list)
        icd_transform = {}
        if len(self.coords) == 1:  # was ncoords > 2; changed b/c seems dangerous
            return icd_transform

        for nid, node in sorted(self.nodes.items()):
            cid_d = node.Cd()
            if cid_d:
                nids_transform[cid_d].append(nid)

        nids_all = np.array(sorted(self.point_ids))
        for cid in sorted(nids_transform.keys()):
            nids = np.array(nids_transform[cid])
            icd_transform[cid] = np.where(np.isin(nids_all, nids))[0]
        return nids_all, nids_transform, icd_transform

    def get_displacement_index_transforms(self):
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Returns
        -------
        icd_transform : dict{int cid : (n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        beta_transforms : dict{in:3x3 float ndarray}
            Dictionary from coordinate id to 3 x 3 transformation
            matrix for that coordinate system.

        Examples
        --------
        # assume GRID 1 has a CD=10
        # assume GRID 2 has a CD=10
        # assume GRID 5 has a CD=50
        >>> model.point_ids
        [1, 2, 5]
        >>> icd_transform, beta_transforms = model.get_displacement_index_transforms()
        >>> icd_transform[10]
        [0, 1]
        >>> beta_transforms[10]
        [1., 0., 0.]
        [0., 0., 1.]
        [0., 1., 0.]

        >>> icd_transform[50]
        [2]
        >>> beta_transforms[50]
        [1., 0., 0.]
        [0., 1., 0.]
        [0., 0., 1.]
        """
        self.deprecated('icd_transform, beta_transforms= model.get_displacement_index_transforms()',
                        'nids_all, nids_transform, icd_transform = model.get_displacement_index()', '1.0')
        nids_transform = defaultdict(list)
        icd_transform = {}
        beta_transforms = {}
        if len(self.coords) == 1:  # was ncoords > 2; changed b/c seems dangerous
            return icd_transform, beta_transforms

        for nid, node in sorted(self.nodes.items()):
            cid_d = node.Cd()
            if cid_d:
                nids_transform[cid_d].append(nid)

        nids_all = np.array(sorted(self.point_ids))
        for cid in sorted(nids_transform.keys()):
            nids = np.array(nids_transform[cid])
            icd_transform[cid] = np.where(np.isin(nids_all, nids))[0]
            beta_transforms[cid] = self.coords[cid].beta()
        return icd_transform, beta_transforms

    def _get_card_name(self, lines):
        """
        Returns the name of the card defined by the provided lines

        Parameters
        ----------
        lines : list[str]
            the lines of the card

        Returns
        -------
        card_name : str
            the name of the card
        """
        card_name = lines[0][:8].rstrip('\t, ').split(',')[0].split('\t')[0].strip('*\t ')
        if len(card_name) == 0:
            return None
        if ' ' in card_name or len(card_name) == 0:
            msg = 'card_name=%r\nline=%r in filename=%r is invalid' \
                  % (card_name, lines[0], self.active_filename)
            print(msg)
            raise CardParseSyntaxError(msg)
        return card_name.upper()

    def increase_card_count(self, card_name, count_num=1):
        """
        Used for testing to check that the number of cards going in is the
        same as each time the model is read verifies proper writing of cards

        Parameters
        ----------
        card_name : str
            the card_name -> 'GRID'
        count_num : int, optional
            the amount to increment by (default=1)

        >>> bdf.read_bdf(bdf_filename)
        >>> bdf.card_count['GRID']
        50
        """
        if card_name in self.card_count:
            self.card_count[card_name] += count_num
        else:
            self.card_count[card_name] = count_num


    def _parse_spc1(self, card_name, cards):
        """adds SPC1s"""
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            constraint_id, dofs, node_ids = get_spc1_constraint(card_obj)

            if constraint_id in self.spc1:
                spc1 = self.spc1[constraint_id]
            else:
                spc1 = SPC1(self)
                self.spc1[constraint_id] = spc1
            #spc1 = self.spc1.setdefault(constraint_id, SPC1(self))
            #spc1.add_card(card_obj, comment=comment)
            spc1.add(constraint_id, dofs, node_ids, comment=comment)
        self.increase_card_count(card_name, len(cards))

    def _parse_mpc(self, card_name, cards):
        """adds MPCs"""
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            constraint_id, constraint = get_mpc_constraint(card_obj)

            if constraint_id in self.mpc:
                mpc = self.mpc[constraint_id]
            else:
                mpc = MPC(self)
                self.mpc[constraint_id] = mpc
            #mpc = self.mpc.setdefault(constraint_id, MPC(self))
            #mpc.add_card(card_obj, comment=comment)
            mpc.add(constraint_id, constraint, comment=comment)
        for constraint_id, constraint in self.mpc.items():
            constraint.build()
        self.increase_card_count(card_name, len(cards))

    def _parse_spcadd(self, card_name, cards):
        """adds SPCADDs"""
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            constraint_id, node_ids = get_spcadd_constraint(card_obj)

            if constraint_id in self.spcadd:
                spcadd = self.spcadd[constraint_id]
            else:
                spcadd = SPCADD(self)
                self.spcadd[constraint_id] = spcadd
            #spcadd.add_card(card_obj, comment=comment)
            spcadd.add(constraint_id, node_ids, comment=comment)
        self.increase_card_count(card_name, len(cards))

    def _parse_mpcadd(self, card_name, cards):
        """adds MPCADDs"""
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            constraint_id, node_ids = get_spcadd_constraint(card_obj)

            if constraint_id in self.mpcadd:
                mpcadd = self.mpcadd[constraint_id]
            else:
                mpcadd = MPCADD(self)
                self.mpcadd[constraint_id] = mpcadd
            #mpcadd.add_card(card_obj, comment=comment)
            mpcadd.add(constraint_id, node_ids, comment=comment)
        self.increase_card_count(card_name, len(cards))

    def _parse_spc(self, card_name, cards):
        """SPC"""
        self._parse_spci(card_name, cards, SPC, self.spc)

    def _parse_spcd(self, card_name, cards):
        """SPCD"""
        self._parse_spci(card_name, cards, SPCD, self.spcd)

    def _parse_spci(self, card_name, cards, obj, slot):
        """SPC, SPCD"""
        #ncards = defaultdict(int)
        data_comments = defaultdict(list)
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            for i in [0, 1]:
                data = get_spc_constraint(card_obj, i)
                constraint_id, node_id = data[:2]
                #self.log.debug('constraint_id=%s node_id=%s dofs=%s enforced=%s' % (
                    #constraint_id, node_id, dofs, enforced_motion))
                if node_id is None:
                    continue
                data_comments[constraint_id].append((data, comment))
                comment = ''

            for constraint_id, data_commentsi in data_comments.items():
                instance = obj(self)
                slot[constraint_id] = instance
                instance.allocate({card_name : len(data_commentsi)})
                for data_comment in data_commentsi:
                    data, comment = data_comment
                    constraint_id, node_id, dofs, enforced_motion = data
                    instance.add(constraint_id, node_id, dofs, enforced_motion, comment=comment)

    def _parse_ctetra(self, card_name, card):
        """adds ctetras"""
        self._parse_solid(
            'CTETRA', card, 7,
            ('CTETRA4', self.ctetra4),
            ('CTETRA10', self.ctetra10),
        )

    def _parse_cpenta(self, card_name, card):
        """adds cpentas"""
        self._parse_solid(
            'CPENTA', card, 9,
            ('CPENTA6', self.cpenta6),
            ('CPENTA15', self.cpenta15),
        )

    def _parse_cpyram(self, card_name, card):
        """adds cpyrams"""
        self._parse_solid(
            'CPENTA', card, 8,
            ('CPENTA6', self.cpenta6),
            ('CPENTA15', self.cpenta15),
        )

    def _parse_chexa(self, card_name, card):
        """adds chexas"""
        self._parse_solid(
            'CHEXA', card, 11,
            ('CHEXA8', self.chexa8),
            ('CHEXA20', self.chexa20),
        )

    @staticmethod
    def _cardlines_to_card_obj(card_lines, card_name):
        """makes a BDFCard object"""
        fields = to_fields(card_lines, card_name)
        card = wipe_empty_fields(fields)
        card_obj = BDFCard(card, has_none=False)
        return card_obj

    def _parse_darea(self, card_name, cards):
        """adds dareas"""
        self._parse_multi(card_name, cards, self.darea, [5])

    def _parse_dphase(self, card_name, cards):
        """adds dphases"""
        self._parse_multi(card_name, cards, self.dphase, [5])

    def _parse_cmass4(self, card_name, cards):
        """adds cmass4"""
        self._parse_multi(card_name, cards, self.cmass4, [5])

    def _parse_cdamp4(self, card_name, cards):
        """adds cdamp4"""
        self._parse_multi(card_name, cards, self.cdamp4, [5])

    def _parse_pvisc(self, card_name, cards):
        """adds pvisc"""
        self._parse_multi(card_name, cards, self.pvisc, [5])

    def _parse_pdamp(self, card_name, cards):
        """adds pdamp"""
        self._parse_multi(card_name, cards, self.pdamp, [3, 5, 7])

    def _parse_pmass(self, card_name, cards):
        """adds pmass"""
        self._parse_multi(card_name, cards, self.pmass, [3, 5, 7])

    def _parse_multi(self, card_name: str, cards, card_cls, icard: list[int]):
        """parses a DAREA, DPHASE, CDAMP4, CMASS4, CVISC, PMASS, PDAMP, ???"""
        datas = []
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)

            data = card_cls.parse(card_obj, icard=0, comment=comment)
            datas.append(data)
            for icardi in icard:
                if card_obj.field(icardi):
                    data = card_cls.parse(card_obj, icard=icardi)
                    datas.append(data)

        ncards = len(datas)
        self.increase_card_count(card_name, ncards)
        self.log.debug('  allocating %r' % card_cls.type)
        card_cls.allocate(self.card_count)
        for datai, comment in datas:
            card_cls.add_card(datai, comment)
        self.log.debug('  building %r; n=%s' % (card_cls.type, card_cls.n))
        card_cls.build()

        #def _prepare_darea(self, card, card_obj, comment=''):
            #"""adds a DAREA"""
            ##def add_darea(self, darea, allow_overwrites=False):
                ##key = (darea.sid, darea.p, darea.c)
                ##if key in self.dareas and not allow_overwrites:
                    ##if not darea == self.dareas[key]:
                        ##assert key not in self.dareas, '\ndarea=\n%s oldDArea=\n%s' % (darea, self.dareas[key])
                ##else:
                    ##assert darea.sid > 0
                    ##self.dareas[key] = darea
                    ##self._type_to_id_map[darea.type].append(key)

            #class_instance = DAREA.add_card(card_obj, comment=comment)
            #self.add_darea(class_instance)
            #if card_obj.field(5):
                #class_instance = DAREA.add_card(card_obj, icard=1, comment=comment)
                #self.add_darea(class_instance)

    def _parse_solid(self, card_name, cards, nsplit, pair1, pair2):
        """
        adds the cards to the object

        Parameters
        ----------
        card_name : str
            the card name
        cards : list[(comment, card_obj), ...]
            an series of comments and cards
        nsplit : int >= 0
            the location to identify for a card split (7 for CTETRA4/CTETRA10)
        pair1 : (card_name, slot)
            card_name : str
               the card_name; (e.g., CTETRA4)
            slot : obj
               the place to put the data (e.g., self.ctetra4)
        pair2 : (card_name, slot)
            card_name : str
               the card_name; (e.g., CTETRA10)
            slot : obj
               the place to put the data (e.g., self.ctetra10)
        """
        cards1 = []
        cards2 = []
        ncards1 = 0
        ncards2 = 0
        for comment, card_lines in cards:
            card_obj = self._cardlines_to_card_obj(card_lines, card_name)
            if card_obj.nfields == nsplit:
                ncards1 += 1
                cards1.append((comment, card_obj))
            else:
                ncards2 += 1
                cards2.append((comment, card_obj))

        if ncards1:
            name1, obj1 = pair1
            self.increase_card_count(name1, ncards1)
            self.log.debug('  allocating %r' % obj1.type)
            obj1.allocate(self.card_count)
            for comment, card_obj in cards1:
                obj1.add_card(card_obj, comment)
            self.log.debug('  building %r; n=%s' % (obj1.type, obj1.n))
            obj1.build()

        if ncards2:
            name2, obj2 = pair2
            self.increase_card_count(name2, ncards2)
            self.log.debug('  allocating %r' % obj2.type)
            obj2.allocate(self.card_count)
            for comment, card_obj in cards2:
                obj2.add_card(card_obj, comment)
            self.log.debug('  building %r; n=%s' % (obj2.type, obj2.n))
            obj2.build()

    def _parse_cards(self, cards, card_count):
        """creates card objects and adds the parsed cards to the deck"""
        #print('card_count = %s' % card_count)

        if isinstance(cards, dict):
            # TODO: many others...
            cards_to_get_lengths_of = {
                'CTETRA' : self._parse_ctetra,
                'CPENTA' : self._parse_cpenta,
                'CPYRAM' : self._parse_cpyram,
                'CHEXA' : self._parse_chexa,
                'SPC1' : self._parse_spc1,
                'SPCADD' : self._parse_spcadd,
                'SPC' : self._parse_spc,
                'SPCD' : self._parse_spcd,

                'MPC' : self._parse_mpc,
                'MPCADD' : self._parse_mpcadd,

                'DAREA' : self._parse_darea,
                'DPHASE' : self._parse_dphase,

                #'PELAS' : self._parse_pelas,
                'PVISC' : self._parse_pvisc,
                'PDAMP' : self._parse_pdamp,
                'CMASS4' : self._parse_cmass4,
                'PMASS' : self._parse_pmass,

                #'CDAMP1' : self._parse_cdamp1,
                #'CDAMP2' : self._parse_cdamp2,
                #'CDAMP3' : self._parse_cdamp3,
                #'CDAMP4' : self._parse_cdamp4,
            }
            # self._is_cards_dict = True
            # this is the loop that hits...
            card_names = sorted(list(cards.keys()))
            for card_name in card_names:
                if card_name in cards_to_get_lengths_of:
                    card = cards[card_name]
                    ncards = len(card)
                    method = cards_to_get_lengths_of[card_name]
                    self.log.info('dynamic vectorized parse of n%s = %s' % (card_name, ncards))
                    method(card_name, card)
                    del cards[card_name]
                    continue

            card_name_to_obj_mapper = self.card_name_to_obj
            for card_name in card_names:
                card = cards[card_name]
                ncards = len(card)
                if self.is_reject(card_name):# and card_name not in :
                    self.log.warning('n%s = %s (rejecting)' % (card_name, ncards))
                    #self.log.info('  rejecting card_name = %s' % card_name)
                    for comment, card_lines in card:
                        self.rejects.append([_format_comment(comment)] + card_lines)
                    self.increase_card_count(card_name, count_num=ncards)
                elif card_name in cards_to_get_lengths_of:
                    #raise RuntimeError('this should not happen because we deleted the cards above')
                    continue
                else:
                    ncards = len(card)
                    self.log.info('n%s = %r' % (card_name, ncards))
                    if card_name not in card_name_to_obj_mapper:
                        self.log.debug('  card_name=%r is not vectorized' % card_name)
                        for comment, card_lines in card:
                            self.add_card(card_lines, card_name, comment=comment,
                                          is_list=False, has_none=False)
                        del cards[card_name]
                        continue

                    obj = card_name_to_obj_mapper[card_name]
                    if obj is None:
                        self.log.debug('card_name=%r is not vectorized, but should be' % card_name)
                        for comment, card_lines in card:
                            self.add_card(card_lines, card_name, comment=comment,
                                          is_list=False, has_none=False)
                        del cards[card_name]
                        continue

                    self.increase_card_count(card_name, ncards)
                    obj.allocate(self.card_count)
                    self.log.debug('  allocating %r' % card_name)
                    for comment, card_lines in card:
                        #print('card_lines', card_lines)
                        fields = to_fields(card_lines, card_name)
                        card = wipe_empty_fields(fields)
                        card_obj = BDFCard(card, has_none=False)
                        obj.add_card(card_obj, comment=comment)
                    obj.build()
                    self.log.debug('  building %r; n=%s' % (obj.type, obj.n))
                    del cards[card_name]
                #if self.is_reject(card_name):
                    #self.log.info('    rejecting card_name = %s' % card_name)
                    #for cardi in card:
                        #self.increase_card_count(card_name)
                        #self.rejects.append([cardi[0]] + cardi[1])
                #else:
                    #for comment, card_lines in card:
                        #print('card_lines', card_lines)
                        #self.add_card(card_lines, card_name, comment=comment,
                                      #is_list=False, has_none=False)
        else:
            # list - this is the one that's used in the non-vectorized case
            raise NotImplementedError('dict...')
            #for card in cards:
                #card_name, comment, card_lines = card
                #if card_name is None:
                    #msg = 'card_name = %r\n' % card_name
                    #msg += 'card_lines = %s' % card_lines
                    #raise RuntimeError(msg)
                #if self.is_reject(card_name):
                    #self.reject_card_lines(card_name, card_lines, comment)

        self.coords.build()
        self.elements.build()
        self.properties.build()
        self.materials.build()

    def _parse_dynamic_syntax(self, key):
        """
        Applies the dynamic syntax for %varName

        Parameters
        ----------
        key : str
            the uppercased key

        Returns
        -------
        value : int/float/str
            the dynamic value defined by dict_of_vars

        .. seealso:: :func: `set_dynamic_syntax`
        """
        key = key.strip()[1:]
        self.log.debug("dynamic key = %r" % key)
        #self.dict_of_vars = {'P5':0.5,'ONEK':1000.}
        if key not in self.dict_of_vars:
            msg = "key=%r not found in keys=%s" % (key, self.dict_of_vars.keys())
            raise KeyError(msg)
        return self.dict_of_vars[key]

    #def _is_case_control_deck(self, line):
        #line_upper = line.upper().strip()
        #if 'CEND' in line.upper():
            #raise SyntaxError('invalid Case Control Deck card...CEND...')
        #if '=' in line_upper or ' ' in line_upper:
            #return True
        #for card in self.uniqueBulkDataCards:
            #lenCard = len(card)
            #if card in line_upper[:lenCard]:
                #return False
        #return True

    def _parse_primary_file_header(self, bdf_filename):
        """
        Extract encoding, nastran_format, and punch from the primary BDF.

        Parameters
        ----------
        bdf_filename : str
            the input filename

        ..code-block :: python

            $ pyNastran: version=NX
            $ pyNastran: encoding=latin-1
            $ pyNastran: punch=True
            $ pyNastran: dumplines=True
            $ pyNastran: nnodes=10
            $ pyNastran: nelements=100
            $ pyNastran: skip_cards=PBEAM,CBEAM
            $ pyNastran: units=in,lb,s

        ..warning :: pyNastran lines must be at the top of the file
        """
        with open(bdf_filename, 'r') as bdf_file:
            check_header = True
            while check_header:
                try:
                    line = bdf_file.readline()
                except Exception:
                    break

                if line.startswith('$'):
                    key, value = _parse_pynastran_header(line)

                    if key:
                        #print('pyNastran key=%s value=%s' % (key, value))
                        if key == 'version':
                            self.nastran_format = value
                        elif key == 'encoding':
                            self._encoding = value
                        elif key == 'punch':
                            self.punch = True if value == 'true' else False
                        elif key in ['nnodes', 'nelements']:
                            pass
                        elif key == 'dumplines':
                            self.dumplines = True if value == 'true' else False
                        elif key == 'skip_cards':
                            cards = {value.strip() for value in value.upper().split(',')}
                            self.cards_to_read = self.cards_to_read - cards
                        elif 'skip ' in key:
                            type_to_skip = key[5:].strip()
                            #values = [int(value) for value in value.upper().split(',')]
                            values = parse_patran_syntax(value)
                            if type_to_skip not in self.object_attributes():
                                raise RuntimeError('%r is an invalid key' % type_to_skip)
                            if type_to_skip not in self.values_to_skip:
                                self.values_to_skip[type_to_skip] = values
                            else:
                                self.values_to_skip[type_to_skip] = np.hstack([
                                    self.values_to_skip[type_to_skip],
                                    values
                                ])
                        #elif key == 'skip_elements'
                        #elif key == 'skip_properties'
                        elif key == 'units':
                            self.units = [value.strip() for value in value.upper().split(',')]
                        else:
                            raise NotImplementedError(key)
                    else:
                        break
                else:
                    break

    def _verify_bdf(self, xref=None):
        """
        Cross reference verification method.
        """
        if xref is None:
            xref = self._xref
        #for key, card in sorted(self.params.items()):
            #card._verify(xref)
        for key, card in sorted(self.nodes.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, card in sorted(self.coords.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, card in sorted(self.elements.items()):
            try:
                card._verify(xref)
            except Exception:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                print(repr(traceback.format_exception(exc_type, exc_value,
                                                      exc_traceback)))
                print(str(card))

                #raise
        for key, card in sorted(self.properties.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, card in sorted(self.materials.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise

        for key, card in sorted(self.dresps.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise

        for key, card in sorted(self.dvcrels.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, card in sorted(self.dvmrels.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, card in sorted(self.dvprels.items()):
            try:
                card._verify(xref)
            except Exception:
                print(str(card))
                raise
        for key, cards in sorted(self.dvgrids.items()):
            for card in cards:
                try:
                    card._verify(xref)
                except Exception:
                    print(str(card))
                    raise


IGNORE_COMMENTS = (
    '$EXECUTIVE CONTROL DECK',
    '$CASE CONTROL DECK',
    'NODES', 'SPOINTS', 'EPOINTS', 'ELEMENTS',
    'PARAMS', 'PROPERTIES', 'ELEMENTS_WITH_PROPERTIES',
    'ELEMENTS_WITH_NO_PROPERTIES (PID=0 and unanalyzed properties)',
    'UNASSOCIATED_PROPERTIES',
    'MATERIALS', 'THERMAL MATERIALS',
    'CONSTRAINTS', 'SPCs', 'MPCs', 'RIGID ELEMENTS',
    'LOADS', 'AERO', 'STATIC AERO', 'AERO CONTROL SURFACES',
    'FLUTTER', 'GUST', 'DYNAMIC', 'OPTIMIZATION',
    'COORDS', 'THERMAL', 'TABLES', 'RANDOM TABLES',
    'SETS', 'CONTACT', 'REJECTS', 'REJECT_LINES',
    'PROPERTIES_MASS', 'MASSES')

def _clean_comment(comment):
    """
    Removes specific pyNastran comment lines so duplicate lines aren't
    created.

    Parameters
    ----------
    comment : str
        the comment to possibly remove

    Returns
    -------
    updated_comment : str
        the comment
    """
    if comment == '':
        pass
    elif comment in IGNORE_COMMENTS:
        comment = ''
    elif 'pynastran' in comment.lower():
        comment = ''

    #if comment:
        #print(comment)
    return comment

def _lines_to_decks(lines, punch):
    """
    Splits the lines into their deck.
    """
    executive_control_lines = []
    case_control_lines = []
    bulk_data_lines = []

    if punch:
        bulk_data_lines = lines
    else:
        flag = 1
        for i, line in enumerate(lines):
            if flag == 1:
                #line = line.upper()
                if line.upper().startswith('CEND'):
                    assert flag == 1
                    flag = 2
                executive_control_lines.append(line.rstrip())
            elif flag == 2:
                uline = line.upper()
                if 'BEGIN' in uline and ('BULK' in uline or 'SUPER' in uline):
                    assert flag == 2
                    flag = 3
                case_control_lines.append(line.rstrip())
            else:
                break
        for line in lines[i:]:
            bulk_data_lines.append(line.rstrip())

        _check_valid_deck(flag)

    del lines
    #for line in bulk_data_lines:
        #print(line)

    # clean comments
    executive_control_lines = [_clean_comment(line) for line in executive_control_lines]
    case_control_lines = [_clean_comment(line) for line in case_control_lines]
    return executive_control_lines, case_control_lines, bulk_data_lines

def _check_valid_deck(flag):
    """Crashes if the flag is set wrong"""
    if flag != 3:
        if flag == 1:
            found = ' - Executive Control Deck\n'
            missing = ' - Case Control Deck\n'
            missing += ' - Bulk Data Deck\n'
        elif flag == 2:
            found = ' - Executive Control Deck\n'
            found += ' - Case Control Deck\n'
            missing = ' - Bulk Data Deck\n'
        else:
            raise RuntimeError('flag=%r is not [1, 2, 3]' % flag)

        msg = 'This is not a valid BDF (a BDF capable of running Nastran).\n\n'
        msg += 'The following sections were found:\n%s\n' % found
        msg += 'The following sections are missing:\n%s\n' % missing
        msg += 'If you do not have an Executive Control Deck or a Case Control Deck:\n'
        msg += '  1.  call read_bdf(...) with `punch=True`\n'
        msg += "  2.  Add '$ pyNastran : punch=True' to the top of the main file\n"
        msg += '  3.  Name your file *.pch\n\n'
        msg += 'You cannot read a deck that has an Executive Control Deck, but\n'
        msg += 'not a Case Control Deck (or vice versa), even if you have a Bulk Data Deck.\n'
        raise RuntimeError(msg)
