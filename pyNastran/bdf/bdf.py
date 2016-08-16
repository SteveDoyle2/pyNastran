# coding: utf-8
# pylint: disable=W0201,R0915,R0912
"""
Main BDF class.  Defines:
  - BDF

see https://docs.plm.automation.siemens.com/tdoc/nxnastran/10/help/#uid:index
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
import sys
import traceback
from codecs import open as codec_open
from collections import defaultdict

from six import string_types, iteritems, itervalues, iterkeys
from six.moves.cPickle import load, dump

import numpy as np

from pyNastran.bdf.utils import _parse_pynastran_header
from pyNastran.utils import object_attributes, print_bad_path
from pyNastran.bdf.utils import (to_fields, get_include_filename,
                                 parse_executive_control_deck, parse_patran_syntax)
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_16 import print_card_16
from pyNastran.bdf.cards.utils import wipe_empty_fields

from pyNastran.utils import _filename
from pyNastran.utils.log import get_logger2

from pyNastran.bdf.write_path import write_include
from pyNastran.bdf.bdf_interface.assign_type import (integer,
                                                     integer_or_string, string)

from pyNastran.bdf.cards.elements.elements import CFAST, CGAP, CRAC2D, CRAC3D, PLOTEL
from pyNastran.bdf.cards.properties.properties import PFAST, PGAP, PRAC2D, PRAC3D, PCONEAX
from pyNastran.bdf.cards.properties.solid import PLSOLID, PSOLID, PIHEX, PCOMPS

from pyNastran.bdf.cards.elements.springs import (CELAS1, CELAS2, CELAS3, CELAS4,)
from pyNastran.bdf.cards.properties.springs import PELAS, PELAST

from pyNastran.bdf.cards.elements.solid import (CTETRA, CPYRAM, CPENTA, CHEXA, CIHEX1)
from pyNastran.bdf.cards.elements.rigid import RBAR, RBAR1, RBE1, RBE2, RBE3, RROD, RSPLINE

from pyNastran.bdf.cards.elements.shell import (CQUAD, CQUAD4, CQUAD8, CQUADR, CQUADX,
                                                CSHEAR, CTRIA3, CTRIA6, CTRIAX,
                                                CTRIAX6, CTRIAR)
from pyNastran.bdf.cards.properties.shell import PSHELL, PCOMP, PCOMPG, PSHEAR, PLPLANE
from pyNastran.bdf.cards.elements.bush import CBUSH, CBUSH1D, CBUSH2D
from pyNastran.bdf.cards.properties.bush import PBUSH, PBUSH1D, PBUSHT
from pyNastran.bdf.cards.elements.damper import (CVISC, CDAMP1, CDAMP2, CDAMP3, CDAMP4,
                                                 CDAMP5)
from pyNastran.bdf.cards.properties.damper import PVISC, PDAMP, PDAMP5, PDAMPT
from pyNastran.bdf.cards.elements.rods import CROD, CONROD, CTUBE
from pyNastran.bdf.cards.elements.bars import CBAR, CBEAM3, CBEND
from pyNastran.bdf.cards.elements.beam import CBEAM
from pyNastran.bdf.cards.properties.rods import PROD, PTUBE
from pyNastran.bdf.cards.properties.bars import PBAR, PBARL  # PBEND
from pyNastran.bdf.cards.properties.beam import PBEAM, PBEAML, PBCOMP
# CMASS5
from pyNastran.bdf.cards.elements.mass import CONM1, CONM2, CMASS1, CMASS2, CMASS3, CMASS4
from pyNastran.bdf.cards.properties.mass import PMASS, NSM
from pyNastran.bdf.cards.aero import (AECOMP, AEFACT, AELINK, AELIST, AEPARM, AESTAT,
                                      AESURF, AESURFS, AERO, AEROS, CSSCHD,
                                      CAERO1, CAERO2, CAERO3, CAERO4, CAERO5,
                                      PAERO1, PAERO2, PAERO3, PAERO4, PAERO5,
                                      MONPNT1,
                                      FLFACT, FLUTTER, GUST, MKAERO1,
                                      MKAERO2, SPLINE1, SPLINE2, SPLINE3, SPLINE4,
                                      SPLINE5, TRIM, DIVERG)
from pyNastran.bdf.cards.constraints import (SPC, SPCADD, SPCAX, SPC1,
                                             MPC, MPCADD, SUPORT1, SUPORT, SESUP,
                                             GMSPC)
from pyNastran.bdf.cards.coordinate_systems import (CORD1R, CORD1C, CORD1S,
                                                    CORD2R, CORD2C, CORD2S, CORD3G,
                                                    GMCORD)
from pyNastran.bdf.cards.dmig import DMIG, DMI, DMIJ, DMIK, DMIJI, DMIG_UACCEL
from pyNastran.bdf.cards.deqatn import DEQATN
from pyNastran.bdf.cards.dynamic import (DELAY, DPHASE, FREQ, FREQ1, FREQ2, FREQ4, TSTEP,
                                         TSTEPNL, NLPARM, NLPCI, TF)
from pyNastran.bdf.cards.loads.loads import (
    LSEQ, SLOAD, DAREA, RANDPS, RFORCE, RFORCE1, SPCD, LOADCYN)
from pyNastran.bdf.cards.loads.dloads import ACSRCE, DLOAD, TLOAD1, TLOAD2, RLOAD1, RLOAD2
from pyNastran.bdf.cards.loads.static_loads import (LOAD, GRAV, ACCEL, ACCEL1, FORCE,
                                                    FORCE1, FORCE2, MOMENT, MOMENT1, MOMENT2,
                                                    PLOAD, PLOAD1, PLOAD2, PLOAD4, PLOADX1,
                                                    GMLOAD)

from pyNastran.bdf.cards.materials import (MAT1, MAT2, MAT3, MAT4, MAT5,
                                           MAT8, MAT9, MAT10, MAT11,
                                           MATHP, CREEP, EQUIV)
# TODO: add MATT3, MATT8, MATT9
from pyNastran.bdf.cards.material_deps import MATT1, MATT2, MATT4, MATT5, MATS1

from pyNastran.bdf.cards.methods import EIGB, EIGC, EIGR, EIGP, EIGRL
from pyNastran.bdf.cards.nodes import GRID, GRDSET, SPOINTs, EPOINTs
from pyNastran.bdf.cards.optimization import (DCONADD, DCONSTR, DESVAR, DDVAL, DOPTPRM, DLINK,
                                              DRESP1, DRESP2, DRESP3,
                                              DVCREL1, DVCREL2,
                                              DVMREL1, DVMREL2,
                                              DVPREL1, DVPREL2)
from pyNastran.bdf.cards.params import PARAM
from pyNastran.bdf.cards.bdf_sets import (ASET, BSET, CSET, QSET, USET,
                                          ASET1, BSET1, CSET1, QSET1, USET1,
                                          SET1, SET3, RADSET,
                                          SEBSET, SECSET, SEQSET, # SEUSET
                                          SEBSET1, SECSET1, SEQSET1, # SEUSET1
                                          SESET, SEQSEP)
from pyNastran.bdf.cards.thermal.loads import QBDY1, QBDY2, QBDY3, QHBDY, TEMP, TEMPD, QVOL
from pyNastran.bdf.cards.thermal.thermal import (CHBDYE, CHBDYG, CHBDYP, PCONV, PCONVM,
                                                 PHBDY, CONV, CONVM, RADM, RADBC)
from pyNastran.bdf.cards.bdf_tables import (TABLED1, TABLED2, TABLED3, TABLED4,
                                            TABLEM1, TABLEM2, TABLEM3, TABLEM4,
                                            TABLES1, TABDMP1, TABLEST, TABRND1, TABRNDG, TIC,
                                            DTABLE)
from pyNastran.bdf.cards.contact import BCRPARA, BCTADD, BCTSET, BSURF, BSURFS, BCTPARA
from pyNastran.bdf.case_control_deck import CaseControlDeck
from pyNastran.bdf.bdf_methods import BDFMethods
from pyNastran.bdf.bdf_interface.get_card import GetMethods
from pyNastran.bdf.bdf_interface.add_card import AddMethods
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.bdf_interface.write_mesh import WriteMesh
from pyNastran.bdf.bdf_interface.cross_reference import XrefMesh
from pyNastran.bdf.errors import CrossReferenceError, DuplicateIDsError, CardParseSyntaxError
from pyNastran.bdf.field_writer_16 import print_field_16



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
        ### TODO: remove all the extra methods

        #keys_to_suppress = []
        #method_names = model.object_methods(keys_to_skip=keys_to_suppress)

        #methods_to_remove = [
            #'process_card', 'read_bdf', 'fill_dmigs', 'disable_cards', 'set_dynamic_syntax',
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
            #'add_SEUSET', 'add_SPLINE', 'add_spoint', 'add_TEMPD', 'add_TF', 'add_TRIM',
            #'add_TSTEP', 'add_TSTEPNL', 'add_USET',

            #'add_card', 'add_card_fields', 'add_card_lines', 'add_cmethod', 'add_constraint',
            #'add_constraint_MPC', 'add_constraint_MPCADD',
            #'add_constraint_SPC', 'add_constraint_SPCADD',
            #'add_convection_property', 'add_coord', 'add_creep_material', 'add_damper',
            #'add_dload', 'add_dload_entry', 'add_element', 'add_hyperelastic_material',
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


class BDF(BDFMethods, GetMethods, AddMethods, WriteMesh, XrefMesh):
    """
    NASTRAN BDF Reader/Writer/Editor class.
    """
    #: required for sphinx bug
    #: http://stackoverflow.com/questions/11208997/autoclass-and-instance-attributes
    __slots__ = ['_is_dynamic_syntax']
    def __init__(self, debug=True, log=None, mode='msc'):
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

        self._relpath = True
        if sys.version_info < (2, 6):
            self._relpath = False

        self.log = get_logger2(log, debug)

        #: list of all read in cards - useful in determining if entire BDF
        #: was read & really useful in debugging
        self.card_count = {}
        #: stores the card_count of cards that have been rejected
        self.reject_count = {}

        #: was an ENDDATA card found
        #self.foundEndData = False

        #: allows the BDF variables to be scoped properly (i think...)
        GetMethods.__init__(self)
        AddMethods.__init__(self)
        BDFMethods.__init__(self)
        WriteMesh.__init__(self)
        XrefMesh.__init__(self)

        #: useful in debugging errors in input
        self.debug = debug

        #: flag that allows for OpenMDAO-style optimization syntax to be used
        self._is_dynamic_syntax = False

        #: lines that were rejected b/c they were for a card that isnt supported
        self.reject_lines = []

        #: cards that were created, but not processed
        self.reject_cards = []

        # self.__init_attributes()

        #: the list of possible cards that will be parsed
        self.cards_to_read = set([
            '/',
            'ECHOON', 'ECHOOFF',
            'PARAM',

            ## nodes
            'GRID', 'GRDSET', 'SPOINT', 'EPOINT',
            #'POINT', 'POINTAX', 'RINGAX', 'GRIDG'

            # mass
            'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

            ## elements
            # springs
            'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', # 'CELAS5',
            # bushings
            'CBUSH', 'CBUSH1D', 'CBUSH2D',
            # dampers
            'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',
            'CFAST',

            'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',
            'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
            'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',
            'CTETRA', 'CPYRAM', 'CPENTA', 'CHEXA',
            'CIHEX1',
            'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
            'CGAP',

            ## rigid_elements
            'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3', 'RROD', 'RSPLINE',

            ## plotels
            'PLOTEL',

            ## properties
            'PMASS',
            'PELAS', 'PGAP', 'PFAST', 'PLPLANE',
            'PBUSH', 'PBUSH1D',
            'PDAMP', 'PDAMP5',
            'PROD', 'PBAR', 'PBARL', 'PBEAM', 'PTUBE', 'PBCOMP', # 'PBEND',
            'PBEAML',  # not fully supported
            # 'PBEAM3',

            'PSHELL', 'PCOMP', 'PCOMPG', 'PSHEAR',
            'PSOLID', 'PLSOLID', 'PVISC', 'PRAC2D', 'PRAC3D',
            'PIHEX', 'PCOMPS',
            # PQUAD4

            ## pdampt
            'PDAMPT',

            ## pelast
            'PELAST',

            ## pbusht
            'PBUSHT',

            ## creep_materials
            'CREEP',

            ## materials
            'MAT1', 'MAT2', 'MAT3', 'MAT8', 'MAT9', 'MAT10', 'MAT11', 'MATHP',

            ## Material dependence - MATT1/MATT2/etc.
            'MATT1', 'MATT2', 'MATT4', 'MATT5',  #'MATT3', 'MATT8', 'MATT9',
            'MATS1', #'MATS3', 'MATS8',
            # 'MATHE'
            #'EQUIV', # testing only, should never be activated...

            ## thermal_materials
            'MAT4', 'MAT5',

            ## spcs/spcadds
            'SPC', 'SPCADD', 'SPC1', 'SPCAX',
            'GMSPC',

            ## mpcs/mpcadds
            'MPC', 'MPCADD',

            ## suport/suport1/se_suport
            'SUPORT', 'SUPORT1', 'SESUP',

            ## loads
            'LOAD', 'LSEQ', 'LOADCYN', 'RANDPS',
            'DLOAD', 'SLOAD', 'ACSRCE', 'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
            'FORCE', 'FORCE1', 'FORCE2',
            'MOMENT', 'MOMENT1', 'MOMENT2',
            'GRAV', 'ACCEL', 'ACCEL1',
            'PLOAD', 'PLOAD1', 'PLOAD2', 'PLOAD4',
            'PLOADX1', 'RFORCE', 'RFORCE1',
            'GMLOAD', 'SPCD',
            #thermal
            'QVOL',

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

            ## coords
            'CORD1R', 'CORD1C', 'CORD1S',
            'CORD2R', 'CORD2C', 'CORD2S',
            'GMCORD',

            # temperature cards
            'TEMP', 'TEMPD',
            'QBDY1', 'QBDY2', 'QBDY3', 'QHBDY',
            'CHBDYE', 'CHBDYG', 'CHBDYP',
            'PCONV', 'PCONVM', 'PHBDY',
            'RADBC', 'CONV',  # 'RADM',

            # ---- dynamic cards ---- #
            'DAREA',  ## dareas
            'DPHASE',  ## dphases
            'DELAY',  ## delays
            'NLPARM',  ## nlparms
            'NLPCI',  ## nlpcis
            'TSTEP',  ## tsteps
            'TSTEPNL',  ## tstepnls
            'TF',  ## transfer_functions

            ## frequencies
            'FREQ', 'FREQ1', 'FREQ2', #'FREQ4',

            # direct matrix input cards
            'DMIG', 'DMIJ', 'DMIJI', 'DMIK', 'DMI',


            # optimization cards
            'DEQATN', 'DTABLE',
            'DCONSTR', 'DESVAR', 'DDVAL', 'DRESP1', 'DRESP2', 'DRESP3',
            'DVCREL1', 'DVCREL2',
            'DVPREL1', 'DVPREL2',
            'DVMREL1', 'DVMREL2',
            'DOPTPRM', 'DLINK', 'DCONADD',
            #'DSCREEN',

            'SET1', 'SET3',  ## sets
            'ASET', 'ASET1',  ## asets
            'BSET', 'BSET1',  ## bsets
            'CSET', 'CSET1',  ## csets
            'QSET', 'QSET1',  ## qsets
            'USET', 'USET1',  ## usets


            # super-element sets
            'SESET',  ## se_sets

            'SEBSET', 'SEBSET1',  ## se_bsets
            'SECSET', 'SECSET1',  ## se_csets
            'SEQSET', 'SEQSET1',  ## se_qsets
            #'SEUSET', 'SEUSET1',  ## se_usets
            'SEQSEP',
            #'RADSET',

            #------------------------------------------------------------------
            ## tables
            #'TABLEHT', 'TABRNDG',
            'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',  # dynamic tables - freq/time loads
            'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',  # material tables - temperature

            # nonlinear elastic temperature dependent materials (e.g. creep)
            # sees TABLES1
            'TABLEST',
            # material tables - stress (MATS1, CREEP, MATHP)
            'TABLES1',


            ## modal damping table - tables_sdamping
            'TABDMP1',

            ## random_tables
            # PSD=func(freq); used by RANDPS card
            'TABRND1',
            # gust for aeroelastic response; used by RANDPS card
            'TABRNDG',

            #------------------------------------------------------------------

            # initial conditions - sid (set ID)
            #'TIC',  (in bdf_tables.py)

            #: methods
            'EIGB', 'EIGR', 'EIGRL',

            #: cMethods
            'EIGC', 'EIGP',

            #: contact
            'BCTPARA',  ## bctpara
            'BCRPARA',  ## bcrpara
            'BCTADD',  ## bctadds
            'BCTSET',  ## bctsets
            'BSURF',  ## bsurf
            'BSURFS',  ## bsurfs

            # other
            'INCLUDE',  # '='
            'ENDDATA',
        ])

        case_control_cards = set(['FREQ', 'GUST', 'MPC', 'SPC', 'NLPARM', 'NSM',
                                  'TEMP', 'TSTEPNL', 'INCLUDE'])
        self._unique_bulk_data_cards = self.cards_to_read.difference(case_control_cards)

        #: / is the delete from restart card
        self.special_cards = ['DEQATN', '/']
        self._make_card_parser()

        if self.is_msc:
            self.set_as_msc()
        elif self.is_nx:
            self.set_as_nx()
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
        del state['_card_parser'], state['_card_parser_b'], state['log']
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
            '_card_parser', '_card_parser_b',
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
        .. note:: loads/spcs (not supported) are tricky because you
                  can't replace cards one-to-one...not sure what to do
        """
        for nid, node in iteritems(replace_model.nodes):
            self.nodes[nid] = node
        for eid, elem in iteritems(replace_model.elements):
            self.elements[eid] = elem
        for eid, elem in iteritems(replace_model.rigid_elements):
            self.rigid_elements[eid] = elem
        for pid, prop in iteritems(replace_model.properties):
            self.properties[pid] = prop
        for mid, mat in iteritems(replace_model.materials):
            self.materials[mid] = mat

        for dvid, desvar in iteritems(replace_model.desvars):
            self.desvars[dvid] = desvar
        for dvid, dvprel in iteritems(replace_model.dvprels):
            self.dvprels[dvid] = dvprel
        for dvid, dvmrel in iteritems(replace_model.dvmrels):
            self.dvmrels[dvid] = dvmrel
        for dvid, dvcrel in iteritems(replace_model.dvcrels):
            self.dvcrels[dvid] = dvcrel

    def disable_cards(self, cards):
        """
        Method for removing broken cards from the reader

        Parameters
        ----------
        cards : List[str]; Set[str]
            a list/set of cards that should not be read

        .. python ::

            bdfModel.disable_cards(['DMIG', 'PCOMP'])
        """
        if cards is None:
            return
        elif isinstance(cards, string_types):
            disable_set = set([cards])
        else:
            disable_set = set(cards)
        self.cards_to_read = self.cards_to_read.difference(disable_set)

    def set_error_storage(self, nparse_errors=100, stop_on_parsing_error=True,
                          nxref_errors=100, stop_on_xref_error=True):
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
        assert isinstance(nparse_errors, int), type(nparse_errors)
        assert isinstance(nxref_errors, int), type(nxref_errors)
        self._nparse_errors = nparse_errors
        self._nxref_errors = nxref_errors
        self._stop_on_parsing_error = stop_on_parsing_error
        self._stop_on_xref_error = stop_on_xref_error

    def validate(self):
        """runs some checks on the input data beyond just type checking"""
        #for eid, elem in sorted(iteritems(model.elements)):
            #elem.validate()
        for nid, node in sorted(iteritems(self.nodes)):
            node.validate()
        for cid, coord in sorted(iteritems(self.coords)):
            coord.validate()
        for eid, elem in sorted(iteritems(self.elements)):
            elem.validate()
        for pid, prop in sorted(iteritems(self.properties)):
            prop.validate()

        for eid, elem in sorted(iteritems(self.rigid_elements)):
            elem.validate()
        for eid, plotel in sorted(iteritems(self.plotels)):
            plotel.validate()
        for eid, mass in sorted(iteritems(self.masses)):
            mass.validate()
        for pid, property_mass in sorted(iteritems(self.properties_mass)):
            property_mass.validate()

        #------------------------------------------------
        for mid, mat in sorted(iteritems(self.materials)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.thermal_materials)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATS1)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATS3)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATS8)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT1)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT2)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT3)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT4)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT5)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT8)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.MATT9)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.creep_materials)):
            mat.validate()
        for mid, mat in sorted(iteritems(self.hyperelastic_materials)):
            mat.validate()

        #------------------------------------------------
        for key, loads in sorted(iteritems(self.loads)):
            for load in loads:
                load.validate()
        for key, tic in sorted(iteritems(self.tics)):
            tic.validate()
        for key, dloads in sorted(iteritems(self.dloads)):
            for dload in dloads:
                dload.validate()
        for key, dload_entries in sorted(iteritems(self.dload_entries)):
            for dload_entry in dload_entries:
                dload_entry.validate()

        #------------------------------------------------
        for key, nlpci in sorted(iteritems(self.nlpcis)):
            nlpci.validate()
        for key, nlparm in sorted(iteritems(self.nlparms)):
            nlparm.validate()
        for key, tstep in sorted(iteritems(self.tsteps)):
            tstep.validate()
        for key, tstepnl in sorted(iteritems(self.tstepnls)):
            tstepnl.validate()
        for key, transfer_functions in sorted(iteritems(self.transfer_functions)):
            for transfer_function in transfer_functions:
                transfer_function.validate()
        for key, delay in sorted(iteritems(self.delays)):
            delay.validate()

        #------------------------------------------------
        if self.aeros is not None:
            self.aeros.validate()
        for caero_id, caero in sorted(iteritems(self.caeros)):
            caero.validate()
        for key, paero in sorted(iteritems(self.paeros)):
            paero.validate()
        for spline_id, spline in sorted(iteritems(self.splines)):
            spline.validate()

        for key, aecomp in sorted(iteritems(self.aecomps)):
            aecomp.validate()
        for key, aefact in sorted(iteritems(self.aefacts)):
            aefact.validate()
        for key, aelinks in sorted(iteritems(self.aelinks)):
            for aelink in aelinks:
                aelink.validate()
        for key, aeparam in sorted(iteritems(self.aeparams)):
            aeparam.validate()
        for key, aesurf in sorted(iteritems(self.aesurf)):
            aesurf.validate()
        for key, aesurfs in sorted(iteritems(self.aesurfs)):
            aesurfs.validate()
        for key, aestat in sorted(iteritems(self.aestats)):
            aestat.validate()
        for key, trim in sorted(iteritems(self.trims)):
            trim.validate()
        for key, diverg in sorted(iteritems(self.divergs)):
            diverg.validate()
        for key, csschd in sorted(iteritems(self.csschds)):
            csschd.validate()
        #self.monitor_points = []

        #------------------------------------------------
        if self.aero is not None:
            self.aero.validate()
        for key, flfact in sorted(iteritems(self.flfacts)):
            flfact.validate()
        for key, flutter in sorted(iteritems(self.flutters)):
            flutter.validate()
        for key, gust in sorted(iteritems(self.gusts)):
            gust.validate()
        #self.mkaeros = []

        #------------------------------------------------
        for key, bcs in sorted(iteritems(self.bcs)):
            for bc in bcs:
                bc.validate()
        for key, phbdy in sorted(iteritems(self.phbdys)):
            phbdy.validate()
        for key, convection_property in sorted(iteritems(self.convection_properties)):
            convection_property.validate()
        for key, tempd in sorted(iteritems(self.tempds)):
            tempd.validate()
        #------------------------------------------------
        for key, bcrpara in sorted(iteritems(self.bcrparas)):
            bcrpara.validate()
        for key, bctadd in sorted(iteritems(self.bctadds)):
            bctadd.validate()
        for key, bctpara in sorted(iteritems(self.bctparas)):
            bctpara.validate()
        for key, bctset in sorted(iteritems(self.bctsets)):
            bctset.validate()
        for key, bsurf in sorted(iteritems(self.bsurf)):
            bsurf.validate()
        for key, bsurfs in sorted(iteritems(self.bsurfs)):
            bsurfs.validate()


        #------------------------------------------------
        for key, suport1 in sorted(iteritems(self.suport1)):
            suport1.validate()
        for suport in self.suport:
            suport.validate()
        for se_suport in self.se_suport:
            se_suport.validate()

        for key, spcs in sorted(iteritems(self.spcs)):
            for spc in spcs:
                spc.validate()
        for key, spcadd in sorted(iteritems(self.spcadds)):
            spcadd.validate()

        for key, mpcs in sorted(iteritems(self.mpcs)):
            for mpc in mpcs:
                mpc.validate()
        for key, mpcadd in sorted(iteritems(self.mpcadds)):
            mpcadd.validate()

        #------------------------------------------------
        for key, darea in sorted(iteritems(self.dareas)):
            darea.validate()
        for key, dphase in sorted(iteritems(self.dphases)):
            dphase.validate()

        for pid, pbusht in sorted(iteritems(self.pbusht)):
            pbusht.validate()
        for pid, pdampt in sorted(iteritems(self.pdampt)):
            pdampt.validate()
        for pid, pelast in sorted(iteritems(self.pelast)):
            pelast.validate()

        for pid, frequency in sorted(iteritems(self.frequencies)):
            frequency.validate()
        #------------------------------------------------
        for key, dmi in sorted(iteritems(self.dmis)):
            dmi.validate()
        for key, dmig in sorted(iteritems(self.dmigs)):
            dmig.validate()
        for key, dmij in sorted(iteritems(self.dmijs)):
            dmij.validate()
        for key, dmiji in sorted(iteritems(self.dmijis)):
            dmiji.validate()
        for key, dmik in sorted(iteritems(self.dmiks)):
            dmik.validate()
        #------------------------------------------------
        for key, sets in sorted(iteritems(self.sets)):
            sets.validate()
        for key, uset in sorted(iteritems(self.usets)):
            for useti in uset:
                useti.validate()

        for aset in self.asets:
            aset.validate()
        for bset in self.bsets:
            bset.validate()
        for cset in self.csets:
            cset.validate()
        for qset in self.qsets:
            qset.validate()

        for key, se_set in sorted(iteritems(self.se_sets)):
            se_set.validate()
        for key, se_uset in sorted(iteritems(self.se_usets)):
            se_uset.validate()
        for se_bset in self.se_bsets:
            se_bset.validate()
        for se_cset in self.se_csets:
            se_cset.validate()
        for se_qset in self.se_qsets:
            se_qset.validate()
        #------------------------------------------------
        for key, table in sorted(iteritems(self.tables)):
            table.validate()
        for key, random_table in sorted(iteritems(self.random_tables)):
            random_table.validate()
        for key, table_sdamping in sorted(iteritems(self.tables_sdamping)):
            table_sdamping.validate()
        #------------------------------------------------
        for key, method in sorted(iteritems(self.methods)):
            method.validate()
        for key, cmethod in sorted(iteritems(self.cMethods)):
            cmethod.validate()
        #------------------------------------------------
        for key, dconadd in sorted(iteritems(self.dconadds)):
            dconadd.validate()
        for key, dconstrs in sorted(iteritems(self.dconstrs)):
            for dconstr in dconstrs:
                dconstr.validate()
        for key, desvar in sorted(iteritems(self.desvars)):
            desvar.validate()
        for key, ddval in sorted(iteritems(self.ddvals)):
            ddval.validate()
        for key, dlink in sorted(iteritems(self.dlinks)):
            dlink.validate()
        for key, dresp in sorted(iteritems(self.dresps)):
            dresp.validate()

        if self.dtable is not None:
            self.dtable.validate()
        if self.doptprm is not None:
            self.doptprm.validate()
        for key, dequation in sorted(iteritems(self.dequations)):
            dequation.validate()
        for key, dvprel in sorted(iteritems(self.dvprels)):
            dvprel.validate()
        for key, dvmrel in sorted(iteritems(self.dvmrels)):
            dvmrel.validate()
        for key, dvcrel in sorted(iteritems(self.dvcrels)):
            dvcrel.validate()
        for key, dscreen in sorted(iteritems(self.dscreen)):
            dscreen.validate()
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

        self.log.debug('---starting BDF.read_bdf of %s---' % self.bdf_filename)
        executive_control_lines, case_control_lines, \
            bulk_data_lines = self._get_lines(self.bdf_filename, self.punch)

        self.case_control_lines = case_control_lines
        self.executive_control_lines = executive_control_lines

        sol, method, isol_line = parse_executive_control_deck(executive_control_lines)
        self.update_solution(sol, method, isol_line)

        self.case_control_deck = CaseControlDeck(self.case_control_lines, self.log)
        self.case_control_deck.solmap_toValue = self._solmap_to_value
        self.case_control_deck.rsolmap_toStr = self.rsolmap_toStr

        if self._is_cards_dict:
            cards, card_count = self.get_bdf_cards_dict(bulk_data_lines)
        else:
            cards, card_count = self.get_bdf_cards(bulk_data_lines)
        self._parse_cards(cards, card_count)

        if self.values_to_skip:
            for key, values in iteritems(self.values_to_skip):
                dict_values = getattr(self, key)
                if not isinstance(dict_values, dict):
                    msg = '%r is an invalid type; only dictionaries are supported' % key
                    raise TypeError(msg)
                for value in values:
                    del dict_values[value]
            # TODO: redo get_card_ids_by_card_types & card_count

        self.pop_parse_errors()
        self.fill_dmigs()

        if validate:
            self.validate()

        self.cross_reference(xref=xref)
        self._xref = xref

        self.log.debug('---finished BDF.read_bdf of %s---' % self.bdf_filename)
        self.pop_xref_errors()

    def _read_bdf_helper(self, bdf_filename, encoding, punch, read_includes):
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

        if not os.path.exists(bdf_filename):
            msg = 'cannot find bdf_filename=%r\n%s' % (bdf_filename, print_bad_path(bdf_filename))
            raise IOError(msg)
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
        for name, card_comments in iteritems(self._dmig_temp):
            card0, comment0 = card_comments[0]
            card_name = card0[0]
            card_name = card_name.rstrip(' *').upper()

            if card_name == 'DMIG':
                # if field2 == 'UACCEL':  # special DMIG card
                card = self.dmigs[name]
            elif card_name == 'DMI':
                card = self.dmis[name]
            elif card_name == 'DMIJ':
                card = self.dmijs[name]
            elif card_name == 'DMIJI':
                card = self.dmijis[name]
            elif card_name == 'DMIK':
                card = self.dmiks[name]
            else:
                raise NotImplementedError(card_name)

            for (card_obj, comment) in card_comments:
                card._add_column(card_obj, comment=comment)
            card.finalize()

        self._dmig_temp = defaultdict(list)

    def pop_parse_errors(self):
        """raises an error if there are parsing errors"""
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
                    msg += 'old_element=\n%s\n' % self.elements[eid].print_repr_card()
                    msg += 'new_elements=\n'
                    for elem, eidi in zip(self._duplicate_elements, duplicate_eids):
                        if eidi == eid:
                            msg += elem.print_repr_card()
                    msg += '\n'
                    is_error = True
                    raise DuplicateIDsError(msg)

            if self._duplicate_properties:
                duplicate_pids = [prop.pid for prop in self._duplicate_properties]
                uduplicate_pids = np.unique(duplicate_pids)
                msg += 'self.properties IDs are not unique=%s\n' % uduplicate_pids
                for pid in duplicate_pids:
                    msg += 'old_property=\n%s\n' % self.properties[pid].print_repr_card()
                    msg += 'new_properties=\n'
                    for prop, pidi in zip(self._duplicate_properties, duplicate_pids):
                        if pidi == pid:
                            msg += prop.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_masses:
                duplicate_eids = [elem.eid for elem in self._duplicate_masses]
                uduplicate_eids = np.unique(duplicate_eids)
                msg += 'self.massses IDs are not unique=%s\n' % uduplicate_eids
                for eid in uduplicate_eids:
                    msg += 'old_mass=\n%s\n' % self.masses[eid].print_repr_card()
                    msg += 'new_masses=\n'
                    for elem, eidi in zip(self._duplicate_masses, duplicate_eids):
                        if eidi == eid:
                            msg += elem.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_material=\n%s\n' % self.materials[mid].print_repr_card()
                    msg += 'new_materials=\n'
                    for mat, midi in zip(self._duplicate_materials, duplicate_mids):
                        if midi == mid:
                            msg += mat.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_thermal_materials:
                duplicate_mids = [mat.mid for mat in self._duplicate_thermal_materials]
                uduplicate_mids = np.unique(duplicate_mids)
                msg += 'self.thermal_materials IDs are not unique=%s\n' % uduplicate_mids
                for mid in uduplicate_mids:
                    msg += 'old_thermal_material=\n%s\n' % (
                        self.thermal_materials[mid].print_repr_card())
                    msg += 'new_thermal_materials=\n'
                    for mat, midi in zip(self._duplicate_thermal_materials, duplicate_mids):
                        if midi == mid:
                            msg += mat.print_repr_card()
                    msg += '\n'
                    is_error = True

            if self._duplicate_coords:
                duplicate_cids = [coord.cid for coord in self._duplicate_coords]
                uduplicate_cids = np.unique(duplicate_cids)
                msg += 'self.coords IDs are not unique=%s\n' % uduplicate_cids
                for cid in uduplicate_cids:
                    msg += 'old_coord=\n%s\n' % self.coords[cid].print_repr_card()
                    msg += 'new_coords=\n'
                    for coord, cidi in zip(self._duplicate_coords, duplicate_cids):
                        if cidi == cid:
                            msg += coord.print_repr_card()
                    msg += '\n'
                    is_error = True

            if is_error:
                msg = 'There are dupliate cards.\n\n' + msg

            if self._stop_on_xref_error:
                msg += 'There are parsing errors.\n\n'
                for (card, an_error) in self._stored_parse_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    msg += 'xref errror: %s\n\n'% an_error[0]
                    is_error = True

            if is_error:
                print('%s' % msg)
                raise DuplicateIDsError(msg.rstrip())

    def pop_xref_errors(self):
        """raises an error if there are cross-reference errors"""
        is_error = False
        if self._stop_on_xref_error:
            if self._ixref_errors == 1 and self._nxref_errors == 0:
                raise
            if self._stored_xref_errors:
                msg = 'There are cross-reference errors.\n\n'
                for (card, an_error) in self._stored_xref_errors:
                    msg += '%scard=%s\n' % (an_error[0], card)
                    is_error = True

                if is_error and self._stop_on_xref_error:
                    raise CrossReferenceError(msg.rstrip())

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
                    cards.append([old_card_name, full_comment, card_lines])

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
                        full_comment += backup_comment + '$' + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += '$' + comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += '$' + comment + '\n'
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
            cards.append([old_card_name, backup_comment + full_comment, card_lines])
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
                        full_comment += backup_comment + '$' + comment + '\n'
                    else:
                        full_comment += backup_comment
                    backup_comment = ''
                elif comment:
                    full_comment += '$' + comment + '\n'
                    backup_comment = ''

            elif comment:
                backup_comment += '$' + comment + '\n'
                #print('add backup=%r' % backup_comment)
            #elif comment:
                #backup_comment += '$' + comment + '\n'

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

    def update_solution(self, sol, method, isol_line):
        """
        Updates the overall solution type (e.g. 101,200,600)

        Parameters
        ----------
        sol : int
            the solution type (101, 103, etc)
        method : str
            the solution method (only for SOL=600)
        isol_line : int
            the line to put the SOL/method on
        """
        self.iSolLine = isol_line
        # the integer of the solution type (e.g. SOL 101)
        if sol is None:
            self.sol = None
            self.solMethod = None
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
            self.solMethod = method.strip()
            self.log.debug("sol=%s method=%s" % (self.sol, self.solMethod))
        else:  # very common
            self.solMethod = None

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

        .. note:: Case sensitivity is supported.
        .. note:: Variables should be 7 characters or less to fit in an
          8-character field.

        .. warning:: Type matters!
        """
        self.dict_of_vars = {}
        assert len(dict_of_vars) > 0, 'nvars = %s' % len(dict_of_vars)
        for (key, value) in sorted(iteritems(dict_of_vars)):
            assert len(key) <= 7, ('max length for key is 7; '
                                   'len(%s)=%s' % (key, len(key)))
            assert len(key) >= 1, ('min length for key is 1; '
                                   'len(%s)=%s' % (key, len(key)))
            if not isinstance(key, string_types):
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

    def process_card(self, card_lines):
        """
        Converts card_lines into a card.
        Considers dynamic syntax and removes empty fields

        Parameters
        ----------
        card_lines : List[str]
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
            >>> fields, card_name = model.process_card(card_lines)
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
        self._increase_card_count(card_name)
        if card_name in ['DEQATN']:
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
        self._increase_card_count(card_name)
        if card_name in ['DEQATN']:
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
        self._increase_card_count(card_name)
        if card_name in ['DEQATN']:
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
        class Crash(object):
            def __init__(self):
                pass
            @classmethod
            def add_card(cls, card, comment=''):
                raise NotImplementedError(card)

        self._card_parser = {
            #'=' : (Crash, None),
            '/' : (Crash, None),
            'GRID' : (GRID, self.add_node),
            'SPOINT' : (SPOINTs, self.add_spoint),
            'EPOINT' : (EPOINTs, self.add_epoint),

            'PARAM' : (PARAM, self.add_PARAM),

            'CORD2R' : (CORD2R, self.add_coord),
            'CORD2C' : (CORD2C, self.add_coord),
            'CORD2S' : (CORD2S, self.add_coord),
            'GMCORD' : (GMCORD, self.add_coord),

            'PLOTEL' : (PLOTEL, self.add_plotel),

            'CONROD' : (CONROD, self.add_element),
            'CROD' : (CROD, self.add_element),
            'PROD' : (PROD, self.add_property),
            'CTUBE' : (CTUBE, self.add_element),
            'PTUBE' : (PTUBE, self.add_property),

            'CBAR' : (CBAR, self.add_element),
            'PBAR' : (PBAR, self.add_property),
            'PBARL' : (PBARL, self.add_property),

            'CBEAM' : (CBEAM, self.add_element),
            'PBEAM' : (PBEAM, self.add_property),
            'PBEAML' : (PBEAML, self.add_property),
            'PBCOMP' : (PBCOMP, self.add_property),

            'CBEAM3' : (CBEAM3, self.add_element),
            #'PBEAM3' : (PBEAM3, self.add_property),

            'CBEND' : (CBEND, self.add_element),
            #'PBEND' : (PBEND, self.add_property),

            'CTRIA3' : (CTRIA3, self.add_element),
            'CQUAD4' : (CQUAD4, self.add_element),
            'CQUAD' : (CQUAD, self.add_element),
            'CQUAD8' : (CQUAD8, self.add_element),
            'CQUADX' : (CQUADX, self.add_element),
            'CQUADR' : (CQUADR, self.add_element),
            'CTRIA6' : (CTRIA6, self.add_element),
            'CTRIAR' : (CTRIAR, self.add_element),
            'CTRIAX' : (CTRIAX, self.add_element),
            'CTRIAX6' : (CTRIAX6, self.add_element),
            'PCOMP' : (PCOMP, self.add_property),
            'PCOMPG' : (PCOMPG, self.add_property),
            'PSHELL' : (PSHELL, self.add_property),
            'PLPLANE' : (PLPLANE, self.add_property),

            'CSHEAR' : (CSHEAR, self.add_element),
            'PSHEAR' : (PSHEAR, self.add_property),

            'CTETRA' : (CTETRA, self.add_element),
            'CPYRAM' : (CPYRAM, self.add_element),
            'CPENTA' : (CPENTA, self.add_element),
            'CHEXA' : (CHEXA, self.add_element),
            'CIHEX1' : (CIHEX1, self.add_element),
            'PIHEX' : (PIHEX, self.add_property),
            'PSOLID' : (PSOLID, self.add_property),
            'PLSOLID' : (PLSOLID, self.add_property),
            'PCOMPS' : (PCOMPS, self.add_property),

            'CELAS1' : (CELAS1, self.add_element),
            'CELAS2' : (CELAS2, self.add_element),
            'CELAS3' : (CELAS3, self.add_element),
            'CELAS4' : (CELAS4, self.add_element),
            'CVISC' : (CVISC, self.add_element),
            'PELAST' : (PELAST, self.add_PELAST),

            'CDAMP1' : (CDAMP1, self.add_damper),
            'CDAMP2' : (CDAMP2, self.add_damper),
            'CDAMP3' : (CDAMP3, self.add_damper),
            # CDAMP4 added later because the documentation is wrong
            'CDAMP5' : (CDAMP5, self.add_damper),
            'PDAMP5' : (PDAMP5, self.add_property),

            'CFAST' : (CFAST, self.add_damper),
            'PFAST' : (PFAST, self.add_property),

            'CGAP' : (CGAP, self.add_element),
            'PGAP' : (PGAP, self.add_property),

            'CBUSH' : (CBUSH, self.add_damper),
            'CBUSH1D' : (CBUSH1D, self.add_damper),
            'CBUSH2D' : (CBUSH2D, self.add_damper),
            'PBUSH' : (PBUSH, self.add_property),
            'PBUSH1D' : (PBUSH1D, self.add_property),

            'CRAC2D' : (CRAC2D, self.add_element),
            'PRAC2D' : (PRAC2D, self.add_property),

            'CRAC3D' : (CRAC3D, self.add_element),
            'PRAC3D' : (PRAC3D, self.add_property),

            'PDAMPT' : (PDAMPT, self.add_PDAMPT),
            'PBUSHT' : (PBUSHT, self.add_PBUSHT),

            'PCONEAX' : (PCONEAX, self.add_property),

            'RBAR' : (RBAR, self.add_rigid_element),
            'RBAR1' : (RBAR1, self.add_rigid_element),
            'RBE1' : (RBE1, self.add_rigid_element),
            'RBE2' : (RBE2, self.add_rigid_element),
            'RBE3' : (RBE3, self.add_rigid_element),
            'RROD' : (RROD, self.add_rigid_element),
            'RSPLINE' : (RSPLINE, self.add_rigid_element),


            ## there is no MAT6 or MAT7
            'MAT1' : (MAT1, self.add_structural_material),
            'MAT2' : (MAT2, self.add_structural_material),
            'MAT3' : (MAT3, self.add_structural_material),
            'MAT8' : (MAT8, self.add_structural_material),
            'MAT9' : (MAT9, self.add_structural_material),
            'MAT10' : (MAT10, self.add_structural_material),
            'MAT11' : (MAT11, self.add_structural_material),
            'EQUIV' : (EQUIV, self.add_structural_material),

            ##'MATHE' : (MATHE, self.add_hyperelastic_material),
            'MATHP' : (MATHP, self.add_hyperelastic_material),
            'MAT4' : (MAT4, self.add_thermal_material),
            'MAT5' : (MAT5, self.add_thermal_material),

            'MATS1' : (MATS1, self.add_material_dependence),
            #'MATS3' : (MATS3, self.add_material_dependence),
            #'MATS8' : (MATS8, self.add_material_dependence),
            'MATT1' : (MATT1, self.add_material_dependence),
            'MATT2' : (MATT2, self.add_material_dependence),
            #'MATT3' : (MATT3, self.add_material_dependence),
            'MATT4' : (MATT4, self.add_material_dependence),
            'MATT5' : (MATT5, self.add_material_dependence),
            #'MATT8' : (MATT8, self.add_material_dependence),
            #'MATT9' : (MATT9, self.add_material_dependence),

            ## hasnt been verified, links up to MAT1, MAT2, MAT9 w/ same MID
            'CREEP' : (CREEP, self.add_creep_material),

            'CONM1' : (CONM1, self.add_mass),
            'CONM2' : (CONM2, self.add_mass),
            'CMASS1' : (CMASS1, self.add_mass),
            'CMASS2' : (CMASS2, self.add_mass),
            'CMASS3' : (CMASS3, self.add_mass),
            ## CMASS4 - added later because documentation is wrong

            'MPC' : (MPC, self.add_constraint_MPC),
            'MPCADD' : (MPCADD, self.add_constraint_MPC),

            'SPC' : (SPC, self.add_constraint_SPC),
            'SPC1' : (SPC1, self.add_constraint_SPC),
            'SPCAX' : (SPCAX, self.add_constraint_SPC),
            'SPCADD' : (SPCADD, self.add_constraint_SPC),
            'GMSPC' : (GMSPC, self.add_constraint_SPC),

            'SESUP' : (SESUP, self.add_sesuport), # pseudo-constraint
            'SUPORT' : (SUPORT, self.add_suport), # pseudo-constraint
            'SUPORT1' : (SUPORT1, self.add_suport1),  # pseudo-constraint

            'FORCE' : (FORCE, self.add_load),
            'FORCE1' : (FORCE1, self.add_load),
            'FORCE2' : (FORCE2, self.add_load),
            'MOMENT' : (MOMENT, self.add_load),
            'MOMENT1' : (MOMENT1, self.add_load),
            'MOMENT2' : (MOMENT2, self.add_load),

            'LSEQ' : (LSEQ, self.add_LSEQ),
            'LOAD' : (LOAD, self.add_load),
            'LOADCYN' : (LOADCYN, self.add_load),
            'GRAV' : (GRAV, self.add_load),
            'ACCEL' : (ACCEL, self.add_load),
            'ACCEL1' : (ACCEL1, self.add_load),
            'PLOAD' : (PLOAD, self.add_load),
            'PLOAD1' : (PLOAD1, self.add_load),
            'PLOAD2' : (PLOAD2, self.add_load),
            'PLOAD4' : (PLOAD4, self.add_load),
            'PLOADX1' : (PLOADX1, self.add_load),
            'RFORCE' : (RFORCE, self.add_load),
            'RFORCE1' : (RFORCE1, self.add_load),
            'SLOAD' : (SLOAD, self.add_load),
            'RANDPS' : (RANDPS, self.add_load),
            'GMLOAD' : (GMLOAD, self.add_load),
            'SPCD' : (SPCD, self.add_load),  # enforced displacement
            'QVOL' : (QVOL, self.add_load),  # thermal

            'DLOAD' : (DLOAD, self.add_dload),
            'ACSRCE' : (ACSRCE, self.add_dload_entry),
            'TLOAD1' : (TLOAD1, self.add_dload_entry),
            'TLOAD2' : (TLOAD2, self.add_dload_entry),
            'RLOAD1' : (RLOAD1, self.add_dload_entry),
            'RLOAD2' : (RLOAD2, self.add_dload_entry),

            'FREQ' : (FREQ, self.add_FREQ),
            'FREQ1' : (FREQ1, self.add_FREQ),
            'FREQ2' : (FREQ2, self.add_FREQ),
            'FREQ4' : (FREQ4, self.add_FREQ),

            'DOPTPRM' : (DOPTPRM, self._add_doptprm),
            'DESVAR' : (DESVAR, self.add_DESVAR),
            # BCTSET

            'TEMP' : (TEMP, self.add_thermal_load),
            'QBDY1' : (QBDY1, self.add_thermal_load),
            'QBDY2' : (QBDY2, self.add_thermal_load),
            'QBDY3' : (QBDY3, self.add_thermal_load),
            'QHBDY' : (QHBDY, self.add_thermal_load),
            'PHBDY' : (PHBDY, self.add_PHBDY),

            'CHBDYE' : (CHBDYE, self.add_thermal_element),
            'CHBDYG' : (CHBDYG, self.add_thermal_element),
            'CHBDYP' : (CHBDYP, self.add_thermal_element),
            'PCONV' : (PCONV, self.add_convection_property),
            'PCONVM' : (PCONVM, self.add_convection_property),

            # aero
            'AECOMP' : (AECOMP, self.add_AECOMP),
            'AEFACT' : (AEFACT, self.add_AEFACT),
            'AELINK' : (AELINK, self.add_AELINK),
            'AELIST' : (AELIST, self.add_AELIST),
            'AEPARM' : (AEPARM, self.add_AEPARM),
            'AESTAT' : (AESTAT, self.add_AESTAT),
            'AESURF' : (AESURF, self.add_AESURF),
            'AESURFS' : (AESURFS, self.add_AESURFS),

            'CAERO1' : (CAERO1, self.add_CAERO),
            'CAERO2' : (CAERO2, self.add_CAERO),
            'CAERO3' : (CAERO3, self.add_CAERO),
            'CAERO4' : (CAERO4, self.add_CAERO),
            'CAERO5' : (CAERO5, self.add_CAERO),

            'PAERO1' : (PAERO1, self.add_PAERO),
            'PAERO2' : (PAERO2, self.add_PAERO),
            'PAERO3' : (PAERO3, self.add_PAERO),
            'PAERO4' : (PAERO4, self.add_PAERO),
            'PAERO5' : (PAERO5, self.add_PAERO),

            'SPLINE1' : (SPLINE1, self.add_SPLINE),
            'SPLINE2' : (SPLINE2, self.add_SPLINE),
            'SPLINE3' : (SPLINE3, self.add_SPLINE),
            'SPLINE4' : (SPLINE4, self.add_SPLINE),
            'SPLINE5' : (SPLINE5, self.add_SPLINE),

            # SOL 144
            'AEROS' : (AEROS, self.add_AEROS),
            'TRIM' : (TRIM, self.add_TRIM),
            'DIVERG' : (DIVERG, self.add_DIVERG),

            # SOL 145
            'AERO' : (AERO, self.add_AERO),
            'FLUTTER' : (FLUTTER, self.add_FLUTTER),
            'FLFACT' : (FLFACT, self.add_FLFACT),
            'MKAERO1' : (MKAERO1, self.add_MKAERO),
            'MKAERO2' : (MKAERO2, self.add_MKAERO),

            'GUST' : (GUST, self.add_gust),
            'CSSCHD' : (CSSCHD, self.add_CSSCHD),
            'MONPNT1' : (MONPNT1, self.add_MONPNT),

            #'NLPARM' : (NLPARM, self.add_NLPARM),
            'NLPCI' : (NLPCI, self.add_NLPCI),
            'TSTEP' : (TSTEP, self.add_TSTEP),
            'TSTEPNL' : (TSTEPNL, self.add_TSTEPNL),

            'TF' : (TF, self.add_TF),
            'DELAY' : (DELAY, self.add_DELAY),

            'DCONADD' : (DCONADD, self.add_DCONSTR),
            'DCONSTR' : (DCONSTR, self.add_DCONSTR),
            'DDVAL' : (DDVAL, self.add_DDVAL),
            'DLINK' : (DLINK, self.add_DLINK),

            'DTABLE' : (DTABLE, self.add_DTABLE),
            'DRESP1' : (DRESP1, self.add_DRESP),
            'DRESP2' : (DRESP2, self.add_DRESP), # deqatn
            'DRESP3' : (DRESP3, self.add_DRESP),
            'DVCREL1' : (DVCREL1, self.add_DVCREL),
            'DVCREL2' : (DVCREL2, self.add_DVCREL),
            'DVPREL1' : (DVPREL1, self.add_DVPREL),
            'DVPREL2' : (DVPREL2, self.add_DVPREL), # deqatn
            'DVMREL1' : (DVMREL1, self.add_DVMREL),
            'DVMREL2' : (DVMREL2, self.add_DVMREL), # deqatn

            'TABLED1' : (TABLED1, self.add_table),
            'TABLED2' : (TABLED2, self.add_table),
            'TABLED3' : (TABLED3, self.add_table),
            'TABLED4' : (TABLED4, self.add_table),
            'TABLEM1' : (TABLEM1, self.add_table),
            'TABLEM2' : (TABLEM2, self.add_table),
            'TABLEM3' : (TABLEM3, self.add_table),
            'TABLEM4' : (TABLEM4, self.add_table),

            'TABLES1' : (TABLES1, self.add_table),
            'TABLEST' : (TABLEST, self.add_table),

            'TABDMP1' : (TABDMP1, self.add_table_sdamping),
            'TABRND1' : (TABRND1, self.add_random_table),
            'TABRNDG' : (TABRNDG, self.add_random_table),

            'EIGB' : (EIGB, self.add_method),
            'EIGR' : (EIGR, self.add_method),
            'EIGRL' : (EIGRL, self.add_method),
            'EIGC' : (EIGC, self.add_cmethod),
            'EIGP' : (EIGP, self.add_cmethod),

            'BCRPARA' : (BCRPARA, self.add_BCRPARA),
            'BCTADD' : (BCTADD, self.add_BCTADD),
            'BCTPARA' : (BCTPARA, self.add_BCTPARA),
            'BSURF' : (BSURF, self.add_BSURF),
            'BSURFS' : (BSURFS, self.add_BSURFS),

            'ASET' : (ASET, self.add_ASET),
            'ASET1' : (ASET1, self.add_ASET),

            'BSET' : (BSET, self.add_BSET),
            'BSET1' : (BSET1, self.add_BSET),

            'CSET' : (CSET, self.add_CSET),
            'CSET1' : (CSET1, self.add_CSET),

            'QSET' : (QSET, self.add_QSET),
            'QSET1' : (QSET1, self.add_QSET),

            'USET' : (USET, self.add_USET),
            'USET1' : (USET1, self.add_USET),

            'SET1' : (SET1, self.add_SET),
            'SET3' : (SET3, self.add_SET),

            'SESET' : (SESET, self.add_SESET),

            'SEBSET' : (SEBSET, self.add_SEBSET),
            'SEBSET1' : (SEBSET1, self.add_SEBSET),

            'SECSET' : (SECSET, self.add_SECSET),
            'SECSET1' : (SECSET1, self.add_SECSET),

            'SEQSET' : (SEQSET, self.add_SEQSET),
            'SEQSET1' : (SEQSET1, self.add_SEQSET),

            #'SESUP' : (SESUP, self.add_SESUP),  # pseudo-constraint

            #'SEUSET' : (SEUSET, self.add_SEUSET),
            #'SEUSET1' : (SEUSET1, self.add_SEUSET),

            'NLPARM' : (NLPARM, self.add_NLPARM),
            # BCTSET
        }
        self._card_parser_prepare = {
            'CORD1R' : self._prepare_cord1r,
            'CORD1C' : self._prepare_cord1c,
            'CORD1S' : self._prepare_cord1s,
            #'CORD3G' : self._prepare_CORD3G,

            'DAREA' : self._prepare_darea,
            'DPHASE' : self._prepare_dphase,
            'PMASS' : self._prepare_pmass,
            'CMASS4' : self._prepare_cmass4,
            'CDAMP4' : self._prepare_cdamp4,

            'DMIG' : self._prepare_dmig,
            'DMI' : self._prepare_dmi,
            'DMIJ' : self._prepare_dmij,
            'DMIK' : self._prepare_dmik,
            'DMIJI' : self._prepare_dmiji,

            'DEQATN' : self._prepare_dequatn,

            'PVISC' : self._prepare_pvisc,
            'PELAS' : self._prepare_pelas,
            'PDAMP' : self._prepare_pdamp,

            'TEMPD' : self._prepare_tempd,
            'CONVM' : self._prepare_convm,
            'CONV' : self._prepare_conv,
            'RADM' : self._prepare_radm,
            'RADBC' : self._prepare_radbc,
            # GRDSET-will be last card to update from _card_parser_prepare
            'GRDSET' : self._prepare_grdset,

            'BCTSET' : self._prepare_bctset,
        }

    def _prepare_bctset(self, card, card_obj, comment=''):
        """adds a GRDSET"""
        card = BCTSET.add_card(card_obj, comment=comment, sol=self.sol)
        self.add_BCTSET(card)

    def _prepare_grdset(self, card, card_obj, comment=''):
        """adds a GRDSET"""
        self.gridSet = GRDSET.add_card(card_obj, comment=comment)

    def _prepare_cdamp4(self, card, card_obj, comment=''):
        """adds a CDAMP4"""
        self.add_damper(CDAMP4.add_card(card_obj, comment=comment))
        if card_obj.field(5):
            self.add_damper(CDAMP4.add_card(card_obj, 1, comment=''))
        return card_obj

    def _prepare_convm(self, card, card_obj, comment=''):
        """adds a CONVM"""
        boundary_condition = CONVM.add_card(card_obj, comment=comment)
        self.add_thermal_BC(boundary_condition, boundary_condition.eid)

    def _prepare_conv(self, card, card_obj, comment=''):
        """adds a CONV"""
        boundary_condition = CONV.add_card(card_obj, comment=comment)
        self.add_thermal_BC(boundary_condition, boundary_condition.eid)

    def _prepare_radm(self, card, card_obj, comment=''):
        """adds a RADM"""
        boundary_condition = RADM.add_card(card, comment=comment)
        self.add_thermal_BC(boundary_condition, boundary_condition.radmid)

    def _prepare_radbc(self, card, card_obj, comment=''):
        """adds a RADBC"""
        boundary_condition = RADBC(card_obj, comment=comment)
        self.add_thermal_BC(boundary_condition, boundary_condition.nodamb)

    def _prepare_tempd(self, card, card_obj, comment=''):
        """adds a TEMPD"""
        self.add_TEMPD(TEMPD.add_card(card_obj, 0, comment=comment))
        if card_obj.field(3):
            self.add_TEMPD(TEMPD.add_card(card_obj, 1, comment=''))
            if card_obj.field(5):
                self.add_TEMPD(TEMPD.add_card(card_obj, 2, comment=''))
                if card_obj.field(7):
                    self.add_TEMPD(TEMPD.add_card(card_obj, 3, comment=''))

    def _add_doptprm(self, doptprm, comment=''):
        """adds a DOPTPRM"""
        self.doptprm = doptprm

    def _prepare_dequatn(self, card, card_obj, comment=''):
        """adds a DEQATN"""
        if hasattr(self, 'test_deqatn') or 1:
            self.add_DEQATN(DEQATN.add_card(card_obj, comment=comment))
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
                self.add_DMIG(card)
            else:
                self._dmig_temp[name].append((card_obj, comment))
        else:
            field2 = integer_or_string(card_obj, 2, 'flag')
            if field2 == 0:
                card = DMIG(card_obj, comment=comment)
                self.add_DMIG(card)
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
        self._prepare_dmix(DMI, self.add_DMI, card_obj, comment=comment)

    def _prepare_dmij(self, card, card_obj, comment=''):
        """adds a DMIJ"""
        self._prepare_dmix(DMIJ, self.add_DMIJ, card_obj, comment=comment)

    def _prepare_dmik(self, card, card_obj, comment=''):
        """adds a DMIK"""
        self._prepare_dmix(DMIK, self.add_DMIK, card_obj, comment=comment)

    def _prepare_dmiji(self, card, card_obj, comment=''):
        """adds a DMIJI"""
        self._prepare_dmix(DMIJI, self.add_DMIJI, card_obj, comment=comment)

    def _prepare_cmass4(self, card, card_obj, comment=''):
        """adds a CMASS4"""
        class_instance = CMASS4.add_card(card_obj, icard=0, comment=comment)
        self.add_mass(class_instance)
        if card_obj.field(5):
            class_instance = CMASS4.add_card(card_obj, icard=1, comment=comment)
            self.add_mass(class_instance)

    def _prepare_pelas(self, card, card_obj, comment=''):
        """adds a PELAS"""
        class_instance = PELAS.add_card(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PELAS.add_card(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)

    def _prepare_pvisc(self, card, card_obj, comment=''):
        """adds a PVISC"""
        class_instance = PVISC.add_card(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PVISC.add_card(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)

    def _prepare_pdamp(self, card, card_obj, comment=''):
        """adds a PDAMP"""
        class_instance = PDAMP.add_card(card_obj, icard=0, comment=comment)
        self.add_property(class_instance)
        if card_obj.field(3):
            class_instance = PDAMP.add_card(card_obj, icard=1, comment=comment)
            self.add_property(class_instance)
        if card_obj.field(5):
            class_instance = PDAMP.add_card(card_obj, icard=2, comment=comment)
            self.add_property(class_instance)
        if card_obj.field(7):
            class_instance = PDAMP.add_card(card_obj, icard=3, comment=comment)
            self.add_property(class_instance)

    def _prepare_pmass(self, card, card_obj, comment=''):
        """adds a PMASS"""
        card_instance = PMASS(card_obj, icard=0, comment=comment)
        self.add_property_mass(card_instance)
        for (i, j) in enumerate([3, 5, 7]):
            if card_obj.field(j):
                card_instance = PMASS(card_obj, icard=i+1, comment=comment)
                self.add_property_mass(card_instance)

    def _prepare_darea(self, card, card_obj, comment=''):
        """adds a DAREA"""
        class_instance = DAREA.add_card(card_obj, comment=comment)
        self.add_DAREA(class_instance)
        if card_obj.field(5):
            class_instance = DAREA.add_card(card_obj, icard=1, comment=comment)
            self.add_DAREA(class_instance)

    def _prepare_dphase(self, card, card_obj, comment=''):
        """adds a DPHASE"""
        class_instance = DPHASE.add_card(card_obj, comment=comment)
        self.add_DPHASE(class_instance)
        #if card_obj.field(5):
            #print('card_obj = ', card_obj)
            #class_instance = DPHASE(card_obj, icard=1, comment=comment)
            #self.add_DPHASE(class_instance)

    def _prepare_cord1r(self, card, card_obj, comment=''):
        """adds a CORD1R"""
        class_instance = CORD1R.add_card(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1R.add_card(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

    def _prepare_cord1c(self, card, card_obj, comment=''):
        """adds a CORD1C"""
        class_instance = CORD1C.add_card(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1C.add_card(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

    def _prepare_cord1s(self, card, card_obj, comment=''):
        """adds a CORD1S"""
        class_instance = CORD1S.add_card(card_obj, comment=comment)
        self.add_coord(class_instance)
        if card_obj.field(5):
            class_instance = CORD1S.add_card(card_obj, icard=1, comment=comment)
            self.add_coord(class_instance)

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

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the card_object is not a card-type object...so not a GRID
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

    def add_card_lines(self, card_lines, card_name, comment='', has_none=True):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_lines: list[str]
            the list of the card fields
            input is list of card_lines -> ['GRID, 1, 2, 3.0, 4.0, 5.0']
        card_name : str
            the card_name -> 'GRID'
        comment : str; default=''
            an optional the comment for the card
        has_none : bool; default=True
            can there be trailing Nones in the card data (e.g. ['GRID, 1, 2, 3.0, 4.0, 5.0, '])
        """
        card_name = card_name.upper()
        card_obj, card = self.create_card_object(card_lines, card_name,
                                                 is_list=False, has_none=has_none)
        self._add_card_helper(card_obj, card, card_name, comment)

    def get_xyz_in_coord(self, cid=0, dtype='float32'):
        """
        Gets the xyz points (including SPOINTS) in the desired coordinate frame

        Parameters
        ----------
        cid : int; default=0
            the desired coordinate system
        dtype : str; default='float32'
            the data type of the xyz coordinates

        Returns
        -------
        xyz : (n, 3) ndarray
            the xyz points in the cid coordinate frame

        .. warning:: doesn't support EPOINTs
        """
        nnodes = len(self.nodes)
        nspoints = 0
        spoints = None
        if self.spoints:
            spoints = self.spoints.points
            nspoints = len(spoints)
        if self.epoints is not None:
            raise NotImplementedError('EPOINTs')

        assert nnodes + nspoints > 0, 'nnodes=%s nspoints=%s' % (nnodes, nspoints)
        xyz_cid0 = np.zeros((nnodes + nspoints, 3), dtype=dtype)
        if cid == 0:
            if nspoints:
                nids = self.nodes.keys()
                newpoints = nids + list(spoints)
                newpoints.sort()
                for i, nid in enumerate(newpoints):
                    if nid not in spoints:
                        node = self.nodes[nid]
                        xyz_cid0[i, :] = node.get_position()
            else:
                for i, (nid, node) in enumerate(sorted(iteritems(self.nodes))):
                    xyz = node.get_position()
                    xyz_cid0[i, :] = xyz
        else:
            if nspoints:
                nids = self.nodes.keys()
                newpoints = nids + list(spoints)
                newpoints.sort()
                for i, nid in enumerate(newpoints):
                    if nid not in spoints:
                        node = self.nodes[nid]
                        xyz_cid0[i, :] = node.get_position_wrt(self, cid)
            else:
                for i, (nid, node) in enumerate(sorted(iteritems(self.nodes))):
                    xyz = node.get_position_wrt(self, cid)
                    xyz_cid0[i, :] = xyz
        return xyz_cid0

    def _add_card_helper(self, card_obj, card, card_name, comment=''):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        card_object : BDFCard()
            the card object representation of card
        card : List[str]
            the fields of the card object; used for rejection and special cards
        card_name : str
            the card_name -> 'GRID'
        comment : str
            an optional the comment for the card
        """
        if self._auto_reject:
            self.reject_cards.append(card)
            print('rejecting processed auto=rejected %s' % card)

        if card_name == 'ECHOON':
            self.echo = True
            return
        elif card_name == 'ECHOOFF':
            self.echo = False
            return

        if self.echo:
            try:
                print(print_card_8(card_obj).rstrip())
            except:
                print(print_card_16(card_obj).rstrip())

        if card_name in self._card_parser:
            card_class, add_card_function = self._card_parser[card_name]
            try:
                class_instance = card_class.add_card(card_obj, comment=comment)
                add_card_function(class_instance)
            except TypeError:
                msg = 'problem adding %s' % card_obj
                raise
                raise TypeError(msg)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                #raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                self._iparse_errors += 1
                #self.log.error(card_obj)
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
                #raise
            #except AssertionError as exception:
                #self.log.error(card_obj)

        elif card_name in self._card_parser_prepare:
            add_card_function = self._card_parser_prepare[card_name]
            try:
                add_card_function(card, card_obj)
            except (SyntaxError, AssertionError, KeyError, ValueError) as exception:
                #raise
                # WARNING: Don't catch RuntimeErrors or a massive memory leak can occur
                #tpl/cc451.bdf
                #raise
                # NameErrors should be caught
                self._iparse_errors += 1
                self.log.error(card_obj)
                var = traceback.format_exception_only(type(exception), exception)
                self._stored_parse_errors.append((card, var))
                if self._iparse_errors > self._nparse_errors:
                    self.pop_parse_errors()
            #except AssertionError as exception:
                #self.log.error(card_obj)
                #raise
        else:
            #raise RuntimeError(card_obj)
            self.reject_cards.append(card_obj)

    def add_card_class(self, class_instance):
        """
        Adds a card object to the BDF object.

        Parameters
        ----------
        class_instance : BaseCard()
            the card class representation of card
        """
        #if card_name == 'ECHOON':
            #self.echo = True
            #return
        #elif card_name == 'ECHOOFF':
            #self.echo = False
            #return

        if self.echo:
            try:
                print(print_card_8(class_instance).rstrip())
            except:
                print(print_card_16(class_instance).rstrip())

        card_name = class_instance.type
        if card_name in self._card_parser:
            add_card_function = self._card_parser[card_name][1]
            add_card_function(class_instance)

        elif card_name in self._card_parser_prepare:
            # TODO: could be faster...
            comment = class_instance.comment
            class_instance.comment = ''
            card_lines = str(class_instance).split('\n')
            self.add_card(card_lines, card_name, comment=comment,
                          is_list=False, has_none=True)
            #add_card_function = self._card_parser_prepare[card_name]
            #add_card_function(card, card_obj)
        else:
            self.reject_cards.append(class_instance)

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

        .. note:: if a card is not supported and not added to the proper
                  lists, this method will fail
        """
        card_stats = [
            'params', 'nodes', 'points', 'elements', 'rigid_elements',
            'properties', 'materials', 'creep_materials',
            'MATT1', 'MATT2', 'MATT3', 'MATT4', 'MATT5', 'MATT8', 'MATT9',
            'MATS1', 'MATS3', 'MATT8',
            'coords', 'mpcs', 'mpcadds',

            # dynamic cards
            'dareas', 'dphases', 'nlparms', 'nlpcis', 'tsteps', 'tstepnls',

            # direct matrix input - DMIG - dict
            'dmis', 'dmigs', 'dmijs', 'dmijis', 'dmiks',
            'dequations',

            # frequencies - dict
            'frequencies',

            # optimization - dict
            'dconadds', 'dconstrs', 'desvars', 'ddvals', 'dlinks', 'dresps',
            'dvcrels', 'dvmrels', 'dvprels',

            # SESETx - dict
            'suport1',
            'se_sets',
            'se_usets',

            # tables
            'tables', 'random_tables',

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
            'gridSet',  # singleton

            'spcs', 'spcadds',

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
            'isStructured', 'uniqueBulkDataCards',
            'nCardLinesMax', 'model_type', 'includeDir',
            'solMethod', 'log',
            'linesPack', 'lineNumbers', 'iSolLine',
            'reject_count', '_relpath', 'isOpened',
            #'foundEndData',
            'specialCards',])

        unsupported_types = ignored_types.union(ignored_types2)
        all_params = object_attributes(self, keys_to_skip=unsupported_types)

        msg = ['---BDF Statistics---']
        # sol
        msg.append('SOL %s\n' % self.sol)

        # loads
        for (lid, loads) in sorted(iteritems(self.loads)):
            msg.append('bdf.loads[%s]' % lid)
            groups_dict = {}
            for load in loads:
                groups_dict[load.type] = groups_dict.get(load.type, 0) + 1
            for name, count_name in sorted(iteritems(groups_dict)):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        # dloads
        for (lid, loads) in sorted(iteritems(self.dloads)):
            msg.append('bdf.dloads[%s]' % lid)
            groups_dict = {}
            for load in loads:
                groups_dict[load.type] = groups_dict.get(load.type, 0) + 1
            for name, count_name in sorted(iteritems(groups_dict)):
                msg.append('  %-8s %s' % (name + ':', count_name))
            msg.append('')

        for (lid, loads) in sorted(iteritems(self.dload_entries)):
            msg.append('bdf.dload_entries[%s]' % lid)
            groups_dict = {}
            for load in loads:
                groups_dict[load.type] = groups_dict.get(load.type, 0) + 1
            for name, count_name in sorted(iteritems(groups_dict)):
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
            groups = set([])

            if not isinstance(card_group, dict):
                msg = '%s is a %s; not dictionary' % (card_group_name, type(card_group))
                raise RuntimeError(msg)
            for card in itervalues(card_group):
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
            for name, counter in sorted(iteritems(self.card_count)):
                if name not in self.cards_to_read:
                    msg.append('  %-8s %s' % (name + ':', counter))
        msg.append('')
        if return_type == 'string':
            return '\n'.join(msg)
        else:
            return msg

    def get_displacement_index_transforms(self):
        """
        Get index and transformation matricies for nodes with
        their output in coordinate systems other than the global.
        Used in combination with ``OP2.transform_displacements_to_global``

        Returns
        ----------
        i_transform : dict{int:(n,) int ndarray}
            Dictionary from coordinate id to index of the nodes in
            ``self.point_ids`` that their output (`CD`) in that
            coordinate system.
        beta_transforms : dict{in:3x3 float ndarray}
            Dictionary from coordinate id to 3 x 3 transformation
            matrix for that coordinate system.

        Example
        -------
        # assume GRID 1 has a CD=10
        # assume GRID 2 has a CD=10
        # assume GRID 5 has a CD=50
        >>> model.point_ids
        [1, 2, 5]
        >>> i_transform, beta_transforms = model.get_displacement_index_transforms()
        >>> i_transform[10]
        [0, 1]
        >>> beta_transforms[10]
        [1., 0., 0.]
        [0., 0., 1.]
        [0., 1., 0.]

        >>> i_transform[50]
        [2]
        >>> beta_transforms[50]
        [1., 0., 0.]
        [0., 1., 0.]
        [0., 0., 1.]
        """
        nids_transform = {}
        i_transform = {}
        beta_transforms = {}
        if len(self.coords) == 1:  # was ncords > 2; changed b/c seems dangerous
            return i_transform, beta_transforms

        for nid, node in sorted(iteritems(self.nodes)):
            cid_d = node.Cd()
            if cid_d:
                if cid_d not in nids_transform:
                    nids_transform[cid_d] = []
                nids_transform[cid_d].append(nid)

        nids_all = np.array(sorted(self.point_ids))
        for cid in sorted(iterkeys(nids_transform)):
            nids = np.array(nids_transform[cid])
            i_transform[cid] = np.where(np.in1d(nids_all, nids))[0]
            beta_transforms[cid] = self.coords[cid].beta()
        return i_transform, beta_transforms

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

    def _show_bad_file(self, bdf_filename):
        """
        Prints the 10 lines before the UnicodeDecodeError occurred.

        Parameters
        ----------
        bdf_filename : str
            the filename to print the lines of
        """
        lines = []
        print('ENCODING - show_bad_file=%r' % self._encoding)

        # rU
        with codec_open(_filename(bdf_filename), 'r', encoding=self._encoding) as bdf_file:
            iline = 0
            nblank = 0
            while 1:
                try:
                    line = bdf_file.readline().rstrip()
                except UnicodeDecodeError:
                    iline0 = max([iline - 10, 0])
                    self.log.error('filename=%s' % self.bdf_filename)
                    for iline1, line in enumerate(lines[iline0:iline]):
                        self.log.error('lines[%i]=%r' % (iline0 + iline1, line))
                    msg = "\n%s encoding error on line=%s of %s; not '%s'" % (
                        self._encoding, iline, bdf_filename, self._encoding)
                    raise RuntimeError(msg)
                if line:
                    nblank = 0
                else:
                    nblank += 1
                if nblank == 20:
                    raise RuntimeError('20 blank lines')
                iline += 1
                lines.append(line)

    def _get_lines(self, bdf_filename, punch=False):
        """
        Opens the bdf and extracts the lines

        Parameters
        ----------

        bdf_filename : str
            the main bdf_filename
        punch : bool, optional
            is this a punch file (default=False; no executive/case control decks)

        Returns
        -------
        executive_control_lines : list[str]
            the executive control deck as a list of strings
        case_control_lines : list[str]
            the case control deck as a list of strings
        bulk_data_lines : list[str]
            the bulk data deck as a list of strings
        """
        #: the directory of the 1st BDF (include BDFs are relative to this one)
        self.include_dir = os.path.dirname(os.path.abspath(bdf_filename))

        with self._open_file(bdf_filename, basename=True) as bdf_file:
            try:
                lines = bdf_file.readlines()
            except:
                self._show_bad_file(bdf_filename)

        nlines = len(lines)

        i = 0
        while i < nlines:
            try:
                line = lines[i].rstrip('\r\n\t')
            except IndexError:
                break
            uline = line.upper()
            if uline.startswith('INCLUDE'):
                j = i + 1
                line_base = line.split('$')[0]
                include_lines = [line_base.strip()]
                # print('----------------------')

                line_base = line_base[8:].strip()
                if line_base.startswith("'") and line_base.endswith("'"):
                    pass
                else:
                    while not line.split('$')[0].endswith("'") and j < nlines:
                        # print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))
                        try:
                            line = lines[j].split('$')[0].strip()
                        except IndexError:
                            # print('bdf_filename=%r' % bdf_filename)
                            crash_name = 'pyNastran_crash.bdf'
                            self._dump_file(crash_name, lines, i+1)
                            msg = 'There was an invalid filename found while parsing (index).\n'
                            msg += 'Check the end of %r\n' % crash_name
                            msg += 'bdf_filename2 = %r' % bdf_filename
                            raise IndexError(msg)
                         #print('endswith_quote=%s; %r' % (
                             #line.split('$')[0].strip().endswith(""), line.strip()))
                        include_lines.append(line.strip())
                        j += 1
                    # print('j=%s nlines=%s less?=%s'  % (j, nlines, j < nlines))

                    #print('*** %s' % line)
                    #bdf_filename2 = line[7:].strip(" '")
                    #include_lines = [line] + lines[i+1:j]
                #print(include_lines)
                bdf_filename2 = get_include_filename(include_lines, include_dir=self.include_dir)
                if self.read_includes:
                    try:
                        self._open_file_checks(bdf_filename2)
                    except IOError:
                        crash_name = 'pyNastran_crash.bdf'
                        self._dump_file(crash_name, lines, j)
                        msg = 'There was an invalid filename found while parsing.\n'
                        msg += 'Check the end of %r\n' % crash_name
                        msg += 'bdf_filename2 = %r\n' % bdf_filename2
                        msg += 'abs_filename2 = %r\n' % os.path.abspath(bdf_filename2)
                        #msg += 'len(bdf_filename2) = %s' % len(bdf_filename2)
                        print(msg)
                        raise
                        raise IOError(msg)

                    with self._open_file(bdf_filename2, basename=False) as bdf_file:
                        #print('bdf_file.name = %s' % bdf_file.name)
                        lines2 = bdf_file.readlines()

                    #print('lines2 = %s' % lines2)
                    nlines += len(lines2)

                    #line2 = lines[j].split('$')
                    #if not line2[0].isalpha():
                        #print('** %s' % line2)

                    include_comment = '\n$ INCLUDE processed:  %s\n' % bdf_filename2
                    #for line in lines2:
                        #print("  ?%s" % line.rstrip())
                    lines = lines[:i] + [include_comment] + lines2 + lines[j:]
                    #for line in lines:
                        #print("  *%s" % line.rstrip())
                else:
                    lines = lines[:i] + lines[j:]
                    self.reject_lines.append(include_lines)
                    #self.reject_lines.append(write_include(bdf_filename2))
            i += 1

        if self.dumplines:
            self._dump_file('pyNastran_dump.bdf', lines, i)

        return _lines_to_decks(lines, i, punch)

    def _dump_file(self, bdf_dump_filename, lines, i):
        """
        Writes a BDF up to some failed line index

        Parameters
        ----------
        bdf_dump_filename : str
            the bdf filename to dump
        lines : List[str]
            the entire list of lines
        i : int
            the last index to write
        """
        with codec_open(_filename(bdf_dump_filename),
                        'w', encoding=self._encoding) as crash_file:
            for line in lines[:i]:
                crash_file.write(line)

    def _increase_card_count(self, card_name, count_num=1):
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

    def _open_file_checks(self, bdf_filename, basename=False):
        """
        Verifies that the BDF about to be opened:
           1.  Exists
           2.  Is Unique
           3.  Isn't an OP2
           4.  Is a File
        """
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)

        if not os.path.exists(_filename(bdf_filename_inc)):
            msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
            msg += 'cwd: %r\n' % os.getcwd()
            msg += 'include_dir: %r\n' % self.include_dir
            msg += print_bad_path(bdf_filename_inc)
            print(msg)
            raise IOError(msg)
        elif bdf_filename_inc.endswith('.op2'):
            print(msg)
            msg = 'Invalid filetype: bdf_filename=%r' % bdf_filename_inc
            raise IOError(msg)
        bdf_filename = bdf_filename_inc

        if bdf_filename in self.active_filenames:
            msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                % (bdf_filename, self.active_filenames)
            print(msg)
            raise RuntimeError(msg)
        elif os.path.isdir(_filename(bdf_filename)):
            current_filename = self.active_filename if len(self.active_filenames) > 0 else 'None'
            msg = 'Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                bdf_filename_inc, current_filename)
            print(msg)
            raise IOError(msg)
        elif not os.path.isfile(_filename(bdf_filename)):
            msg = 'Not a file: bdf_filename=%r' % bdf_filename
            print(msg)
            raise IOError(msg)

    def _open_file(self, bdf_filename, basename=False, check=True):
        """
        Opens a new bdf_filename with the proper encoding and include directory

        Parameters
        ----------
        bdf_filename : str
            the filename to open
        basename : bool (default=False)
            should the basename of bdf_filename be appended to the include directory
        """
        if basename:
            bdf_filename_inc = os.path.join(self.include_dir, os.path.basename(bdf_filename))
        else:
            bdf_filename_inc = os.path.join(self.include_dir, bdf_filename)

        self._validate_open_file(bdf_filename, bdf_filename_inc, check)


        self.log.debug('opening %r' % bdf_filename_inc)
        self.active_filenames.append(bdf_filename_inc)

        #print('ENCODING - _open_file=%r' % self._encoding)
        bdf_file = codec_open(_filename(bdf_filename_inc), 'r', encoding=self._encoding)
        return bdf_file

    def _validate_open_file(self, bdf_filename, bdf_filename_inc, check):
        """
        checks that the file doesn't have obvious errors
         - hasn't been used
         - not a directory
         - is a file

        Parameters
        ----------
        bdf_filename : str
           the current bdf filename
        bdf_filename_inc : str
           the next bdf filename

        Raises
        ------
        RuntimeError : file is active
        IOError : Invalid file type
        """
        if check:
            if not os.path.exists(_filename(bdf_filename_inc)):
                msg = 'No such bdf_filename: %r\n' % bdf_filename_inc
                msg += 'cwd: %r\n' % os.getcwd()
                msg += 'include_dir: %r\n' % self.include_dir
                msg += print_bad_path(bdf_filename_inc)
                raise IOError(msg)
            elif bdf_filename_inc.endswith('.op2'):
                raise IOError('Invalid filetype: bdf_filename=%r' % bdf_filename_inc)

            bdf_filename = bdf_filename_inc
            if bdf_filename in self.active_filenames:
                msg = 'bdf_filename=%s is already active.\nactive_filenames=%s' \
                    % (bdf_filename, self.active_filenames)
                raise RuntimeError(msg)
            elif os.path.isdir(_filename(bdf_filename)):
                current_filename = self.active_filename if len(self.active_filenames) > 0 else 'None'
                raise IOError('Found a directory: bdf_filename=%r\ncurrent_file=%s' % (
                    bdf_filename_inc, current_filename))
            elif not os.path.isfile(_filename(bdf_filename)):
                raise IOError('Not a file: bdf_filename=%r' % bdf_filename)

    def _parse_cards(self, cards, card_count):
        """creates card objects and adds the parsed cards to the deck"""
        #print('card_count = %s' % card_count)
        if isinstance(cards, dict): # self._is_cards_dict = True
            for card_name, card in sorted(iteritems(cards)):
                if self.is_reject(card_name):
                    self.log.info('    rejecting card_name = %s' % card_name)
                    for cardi in card:
                        self._increase_card_count(card_name)
                        self.rejects.append([cardi[0]] + cardi[1])
                else:
                    for comment, card_lines in card:
                        self.add_card(card_lines, card_name, comment=comment,
                                      is_list=False, has_none=False)
        else:
            for card in cards:
                card_name, comment, card_lines = card
                if card_name is None:
                    msg = 'card_name = %r\n' % card_name
                    msg += 'card_lines = %s' % card_lines
                    raise RuntimeError(msg)
                if self.is_reject(card_name):
                    if card_name not in self.card_count:
                        if ' ' in card_name:
                            msg = (
                                'No spaces allowed in card name %r.  '
                                'Should this be a comment?\n%s%s' % (
                                    card_name, comment, card_lines))
                            raise RuntimeError(msg)
                        if card_name in ['SUBCASE ', 'CEND']:
                            raise RuntimeError('No executive/case control deck was defined.')
                        self.log.info('    rejecting card_name = %s' % card_name)
                    self._increase_card_count(card_name)
                    self.rejects.append([comment] + card_lines)
                else:
                    self.add_card(card_lines, card_name, comment=comment,
                                  is_list=False, has_none=False)

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
                except:
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
        #for key, card in sorted(iteritems(self.params)):
            #card._verify(xref)
        for key, card in sorted(iteritems(self.nodes)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.coords)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.elements)):
            try:
                card._verify(xref)
            except Exception as e:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                print(repr(traceback.format_exception(exc_type, exc_value,
                                                      exc_traceback)))
                print(str(card))

                #raise
        for key, card in sorted(iteritems(self.properties)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.materials)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise

        for key, card in sorted(iteritems(self.dresps)):
            #try:
            card._verify(xref)
            #except:
                #print(str(card))
                #raise

        for key, card in sorted(iteritems(self.dvcrels)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.dvmrels)):
            try:
                card._verify(xref)
            except:
                print(str(card))
                raise
        for key, card in sorted(iteritems(self.dvprels)):
            try:
                card._verify(xref)
            except:
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

def _lines_to_decks(lines, i, punch):
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
    del lines
    #for line in bulk_data_lines:
        #print(line)

    # clean comments
    executive_control_lines = [_clean_comment(line) for line in executive_control_lines]
    case_control_lines = [_clean_comment(line) for line in case_control_lines]
    return executive_control_lines, case_control_lines, bulk_data_lines


def main():
    """
    shows off how unicode works becausee it's overly complicated
    """
    import pyNastran
    pkg_path = pyNastran.__path__[0]
    bdf_filename = os.path.abspath(os.path.join(
        pkg_path, '..', 'models', 'solid_bending', 'solid_bending.bdf'))
    model = BDF()
    model.read_bdf(bdf_filename, encoding='latin-1')
    node1 = model.nodes[1]

    # decode when we receive, encode on send
    note = 'helló wörld from two'  # must be same encoding as the header (utf-8)
    #note = b'helló wörld from two\n'.decode('utf-8')
    print(note)

    # this will be wrong because it's inconsistent with the header (utf-8)
    #note = b'á'.decode('latin-1')

    # this will work
    note = 'á'

    # so will this
    #note = b'á'.decode('utf-8')

    # The encoding that goes into the comment must be consisent with the local
    # file, so if your print doesn't work right, your comment will be bad too.
    #
    # If the print is correct, assuming all the characters are supported in your
    # desired encoding, it *should* work.
    print(note)

    # Comments are unmodified, so you can inadvertantly add cards/bugs.
    # A comment is a single string where all lines start with $ and end
    # with an endline character.
    node1.comment = '$ ' + note + '\n'

    # in other words, msg is a bad comment:
    msg = '$ line 1\n'
    msg += 'line 2\n'
    msg += '$ line 3\n'
    print(msg)
    model.write_bdf('test.bdf')


if __name__ == '__main__':  # pragma: no cover
    from pyNastran.bdf.test.test_bdf import main
    #main()
