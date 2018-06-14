from pyNastran.op2.tables.geom.geom1 import GEOM1
from pyNastran.op2.tables.geom.geom2 import GEOM2
from pyNastran.op2.tables.geom.geom3 import GEOM3
from pyNastran.op2.tables.geom.geom4 import GEOM4

from pyNastran.op2.tables.geom.ept import EPT
from pyNastran.op2.tables.geom.mpt import MPT

from pyNastran.op2.tables.geom.dit import DIT
from pyNastran.op2.tables.geom.dynamics import DYNAMICS

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2, FatalError, SortCodeError, DeviceCodeError, FortranMarkerError


def read_op2_geom(op2_filename=None, combine=True, subcases=None,
                  exclude_results=None, include_results=None,
                  validate=True, xref=True,
                  build_dataframe=False, skip_undefined_matrices=True,
                  mode='msc', log=None, debug=True, debug_file=None, encoding=None):
    """
    Creates the OP2 object without calling the OP2 class.

    Parameters
    ----------
    op2_filename : str (default=None -> popup)
        the op2_filename
    combine : bool; default=True
        True : objects are isubcase based
        False : objects are (isubcase, subtitle) based;
                will be used for superelements regardless of the option
    subcases : List[int, ...] / int; default=None->all subcases
        list of [subcase1_ID,subcase2_ID]
    exclude_results / include_results : List[str] / str; default=None
        a list of result types to exclude/include
        one of these must be None
    validate : bool
        runs various checks on the BDF (default=True)
    xref :  bool
        should the bdf be cross referenced (default=True)
    build_dataframe : bool; default=False
        builds a pandas DataFrame for op2 objects
    skip_undefined_matrices : bool; default=False
         True : prevents matrix reading crashes
    log : Log()
        a logging object to write debug messages to
     (.. seealso:: import logging)
    debug : bool; default=False
        enables the debug log and sets the debug in the logger
    debug_file : str; default=None (No debug)
        sets the filename that will be written to
    encoding : str
        the unicode encoding (default=None; system default)

    Returns
    -------
    model : OP2()
        an OP2 object

    .. todo:: creates the OP2 object without all the read methods

    .. note :: this method will change in order to return an object that
               does not have so many methods

    """
    model = OP2Geom(log=log, debug=debug, debug_file=debug_file, mode=mode)
    model.set_subcases(subcases)
    if exclude_results and include_results:
        msg = 'exclude_results or include_results must be None\n'
        msg += 'exclude_results=%r\n' % exclude_results
        msg += 'include_results=%r\n' % include_results
        raise RuntimeError(msg)
    elif exclude_results:
        model.remove_results(exclude_results)
    elif include_results:
        model.set_results(include_results)

    model.read_op2(op2_filename=op2_filename, build_dataframe=build_dataframe,
                   skip_undefined_matrices=skip_undefined_matrices, combine=combine,
                   encoding=encoding)
    if validate:
        model.validate()
    if xref:
        model.cross_reference()
    return model


class OP2GeomCommon(OP2, GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, DIT, DYNAMICS):
    """interface for the OP2Geom class for to loading subclasses"""
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None, mode='msc'):
        """
        Initializes the OP2 object

        Parameters
        ----------
        make_geom : bool; default=False
            reads the BDF tables
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log: log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : default=None -> no debug
            sets the filename that will be written to
        mode : str; default='msc'
            {msc, nx}

        """
        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)
        DIT.__init__(self)
        DYNAMICS.__init__(self)

        OP2.__init__(self, debug, log=log, debug_file=debug_file, mode=mode)
        self.make_geom = True
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altmdtku4.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altd200x7.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mdtku1.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mcso14.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_ds105.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_altcc574.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_adjoint.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_mcso18.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_cqr4optstdis.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_d200ce12.op2
        #: Optimization Table (I think this is NX-specifc)
        self._edom_map = {
            # are these 3 really EDOM?
            #MAT1DOM(103,1,9944)
            #MAT10DOM(2801,28,9945)
            #MODTRAK(6006,60,477)
            (103, 1, 9944) : ['???', self._read_fake],
            (304, 3, 276) : ['???', self._read_fake],
            (404, 4, 277) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],

            (5106,51,471) : ['DCONADD', self._read_fake],
            (4106,41,362) : ['DCONSTR', self._read_fake],
            #DDVAL(7000,70,563)
            #DRESP3(6700,67,433)

            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],
            (504, 5, 246) : ['???', self._read_fake],

            (3106, 31, 352) : ['DESVAR', self._read_fake],
            (3206, 32, 353) : ['DLINK', self._read_fake],
            (3306, 33, 354) : ['DVPREL1', self._read_fake],
            (3406, 34, 355) : ['DVPREL2', self._read_fake],
            #DOPTPRM(4306,43,364)
            (3706, 37, 358) : ['DTABLE', self._read_fake],
            (3806, 38, 359) : ['DRESP1', self._read_fake],
            (3906, 39, 360) : ['DRESP2', self._read_fake],
            (4106, 41, 362) : ['???', self._read_fake],
            (4206, 42, 363) : ['DSCREEN', self._read_fake],
            (4306, 43, 364) : ['???', self._read_fake],
            (4406, 44, 372) : ['DVGRID', self._read_fake],
            #DVSHAP(5006,50,470)
            (5106, 51, 471) : ['???', self._read_fake],
            #DVBSHAP(5806,58,474)
            #DVGEOM(5906,59,356)
            (6006, 60, 477) : ['???', self._read_fake],
            #DRESP3(6700,67,433)
            (6100, 61, 429) : ['DVCREL1', self._read_fake],
            (6200, 62, 430) : ['DVCREL2', self._read_fake],
            (6300, 63, 431) : ['DVMREL1', self._read_fake],
            (6400, 64, 432) : ['DVMREL2', self._read_fake],
            (6006, 60, 477) : ['???', self._read_fake],
            (7000, 70, 563) : ['DCONSTR/DDVAL?', self._read_fake],
        }

        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltsold11b.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conedg01b.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conprop06.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_glueac103a.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_conedg01s.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_sline5.op2
        self._contact_map = {
            (224, 2, 436) : ['???', self._read_fake],
            (724, 7, 441) : ['???', self._read_fake],
            (1224, 12, 446) : ['???', self._read_fake],
            (7110, 71, 588) : ['???', self._read_fake],
            (7210, 72, 589) : ['???', self._read_fake],
            (7410, 74, 591) : ['???', self._read_fake],
            (7510, 75, 592) : ['???', self._read_fake],
            (7710, 77, 594) : ['???', self._read_fake],
            (8110, 81, 598) : ['???', self._read_fake],
            (8301, 83, 605) : ['???', self._read_fake],
            (8710, 87, 449) : ['???', self._read_fake],
            (8810, 88, 603) : ['???', self._read_fake],
            (8920, 89, 614) : ['???', self._read_fake],
            (124, 1, 435) : ['???', self._read_fake],
            (424, 4, 438) : ['???', self._read_fake],
        }
        # F:\Program Files\Siemens\NXNastran\nxn10p1\nxn10p1\nast\tpl\fsw_eng.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_boltld04i.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_eliter17.op2
        # F:\work\pyNastran\pyNastran\master2\pyNastran\bdf\test\nx_spike\out_weld01i.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ac10901a_new.op2
        self._edt_map = {
            (5201, 52, 373) : ['ACMODL', self._read_fake],
            (6301, 63, 397) : ['ADAPT', self._read_fake],
            (7801, 78, 582) : ['AECOMP', self._read_fake],
            (7901, 79, 583) : ['AECOMPL', self._read_fake],
            (7301, 73, 574) : ['AEDW', self._read_fake],
            (4002, 40, 273) : ['AEFACT', self._read_fake],
            (7501, 75, 576) : ['AEFORCE', self._read_fake],
            (2602, 26, 386) : ['AELINK', self._read_fake],
            (2302, 23, 341) : ['AELIST', self._read_fake],
            (7001, 70, 571) : ['AEPARM', self._read_fake],
            (7401, 74, 575) : ['AEPRESS', self._read_fake],
            (3202, 32, 265) : ['AERO', self._read_fake],
            (2202, 22, 340) : ['AEROS', self._read_fake],
            (2102, 21, 339) : ['AESTAT', self._read_fake],
            (2002, 20, 338) : ['AESURF', self._read_fake],
            (7701, 77, 581) : ['AESURFS', self._read_fake],
            (3002, 30, 263) : ['CAERO1', self._read_fake],
            (4301, 43, 167) : ['CAERO2', self._read_fake],
            (4401, 44, 168) : ['CAERO3', self._read_fake],
            (4501, 45, 169) : ['CAERO4', self._read_fake],
            (5001, 50, 175) : ['CAERO5', self._read_fake],
            (6201, 62, 143) : ['CLOAD', self._read_fake],
            (6401, 64, 307) : ['CSSCHD', self._read_fake],
            (104, 1, 81) : ['DEFORM', self._read_fake],
            (2702, 27, 387) : ['DIVERG', self._read_fake],
            (4102, 41, 274) : ['FLFACT', self._read_fake],
            (3902, 39, 272) : ['FLUTTER', self._read_fake],
            (17400, 174, 616) : ['GROUP', self._read_fake],
            (3802, 38, 271) : ['MKAERO1', self._read_fake],
            (3702, 37, 270) : ['MKAERO2', self._read_fake],
            (7601, 76, 577) : ['MONPNT1', self._read_fake],
            (3102, 31, 264) : ['PAERO1', self._read_fake],
            (4601, 46, 170) : ['PAERO2', self._read_fake],
            (4701, 47, 171) : ['PAERO3', self._read_fake],
            (4801, 48, 172) : ['PAERO4', self._read_fake],
            (5101, 51, 176) : ['PAERO5', self._read_fake],
            (5301, 53, 378) : ['PANEL', self._read_fake],
            (3502, 35, 268) : ['SET1', self._read_fake],
            (3602, 36, 269) : ['SET2', self._read_fake],
            (4302, 43, 607) : ['SET3', self._read_fake],
            (3302, 33, 266) : ['SPLINE1', self._read_fake],
            (3402, 34, 267) : ['SPLINE2', self._read_fake],
            (4901, 49, 173) : ['SPLINE3', self._read_fake],
            (6501, 65, 308) : ['SPLINE4', self._read_fake],
            (6601, 66, 309) : ['SPLINE5', self._read_fake],
            (2402, 24, 342) : ['TRIM', self._read_fake],
            (7201, 72, 573) : ['UXVEC', self._read_fake],
            (7108, 822,51) : ['BOLT', self._read_fake],
            (7108, 71, 251) : ['???', self._read_fake],
            (5808, 58, 220) : ['ITER', self._read_fake],
            (14000, 140, 568) : ['SWLDPRM', self._read_fake],
            (11001, 110, 581) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
        }

        # F:\work\pyNastran\examples\Dropbox\move_tpl\beamp10.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ifsr22r.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ifssh22.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ifsr22r.op2
        # F:\work\pyNastran\examples\Dropbox\move_tpl\ifsv02pp.op2
        self._viewtb_map = {
            (10300, 103, 16) : ['QUADP', self._read_fake],
            (10400, 104, 15) : ['TRIAP', self._read_fake],
            (10500, 105, 14) : ['BEAMP', self._read_fake],
            (14100, 141, 18) : ['HEXAP', self._read_fake],
            (14200, 142, 16) : ['PENTAP', self._read_fake],
            (14300, 143, 14) : ['TETRAP', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],
            #(10500, 105, 14) : ['???', self._read_fake],

        }

    def save(self, obj_filename='model.obj', unxref=True):
        # type: (str, bool) -> None
        """Saves a pickleable object"""
        #del self.log
        #del self._card_parser, self._card_parser_prepare

        #print(object_attributes(self, mode="all", keys_to_skip=[]))
        with open(obj_filename, 'wb') as obj_file:
            dump(self, obj_file)

    def _get_table_mapper(self):
        table_mapper = OP2._get_table_mapper(self)

        table_mapper[b'CONTACT'] = [self._read_contact_4, self._read_contact_4]
        table_mapper[b'CONTACTS'] = [self._read_contact_4, self._read_contact_4]
        table_mapper[b'VIEWTB'] = [self._read_viewtb_4, self._read_viewtb_4]
        table_mapper[b'EDT'] = [self._read_edt_4, self._read_edt_4]
        table_mapper[b'EDTS'] = [self._read_edt_4, self._read_edt_4]

        # geometry
        table_mapper[b'GEOM1'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4'] = [self._read_geom4_4, self._read_geom4_4]

        # superelements
        table_mapper[b'GEOM1S'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2S'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3S'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4S'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper[b'GEOM1N'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2N'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3N'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4N'] = [self._read_geom4_4, self._read_geom4_4]

        table_mapper[b'GEOM1OLD'] = [self._read_geom1_4, self._read_geom1_4]
        table_mapper[b'GEOM2OLD'] = [self._read_geom2_4, self._read_geom2_4]
        table_mapper[b'GEOM3OLD'] = [self._read_geom3_4, self._read_geom3_4]
        table_mapper[b'GEOM4OLD'] = [self._read_geom4_4, self._read_geom4_4]
        table_mapper[b'EDOM'] = [self._read_edom4_4, self._read_edom4_4]

        table_mapper[b'EPT'] = [self._read_ept_4, self._read_ept_4]
        table_mapper[b'EPTS'] = [self._read_ept_4, self._read_ept_4]
        table_mapper[b'EPTOLD'] = [self._read_ept_4, self._read_ept_4]

        table_mapper[b'MPT'] = [self._read_mpt_4, self._read_mpt_4]
        table_mapper[b'MPTS'] = [self._read_mpt_4, self._read_mpt_4]

        table_mapper[b'DYNAMIC'] = [self._read_dynamics_4, self._read_dynamics_4]
        table_mapper[b'DYNAMICS'] = [self._read_dynamics_4, self._read_dynamics_4]

        # table objects (e.g. TABLED1)
        table_mapper[b'DIT'] = [self._read_dit_4, self._read_dit_4]
        table_mapper[b'DITS'] = [self._read_dit_4, self._read_dit_4]
        return table_mapper

    def _read_edom4_4(self, data, ndata):
        """
        reads the EDOM table
        SOL 200 design optimization and sensitivity analysis bulk entries.

        """
        return self._read_geom_4(self._edom_map, data, ndata)

    def _read_contact_4(self, data, ndata):
        """
        reads the CONTACT/CONTACTS table
        Table of Bulk Data entry related to surface contact

        """
        return self._read_geom_4(self._contact_map, data, ndata)

    def _read_edt_4(self, data, ndata):
        """
        3.21 EDT
        Aero and element deformations.

        """
        return self._read_geom_4(self._edt_map, data, ndata)

    def _read_viewtb_4(self, data, ndata):
        """
        View information table
        Contains the relationship between each p-element and its view-elements and view-grids.
        """
        return self._read_geom_4(self._viewtb_map, data, ndata)


class OP2Geom(BDF, OP2GeomCommon):
    """creates an interface for the OP2 and BDF classes"""
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None, mode='msc'):
        """
        Initializes the OP2 object

        Parameters
        ----------
        make_geom : bool; default=False
            reads the BDF tables
        debug : bool; default=False
            enables the debug log and sets the debug in the logger
        log: log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : default=None -> no debug
            sets the filename that will be written to
        mode : str; default='msc'
            {msc, nx}

        """
        BDF.__init__(self, debug=debug, log=log)
        OP2GeomCommon.__init__(self, make_geom=make_geom,
                               debug=debug, log=log, debug_file=debug_file, mode=mode)

    def __getstate__(self):
        """clears out a few variables in order to pickle the object"""
        raise NotImplementedError()
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        #adfasd
        state = BDF.__getstate__(self)
        #print(state)
        #state = self.__dict__.copy()

        # Remove the unpicklable entries.
        i = 0
        for key, value in sorted(state.items()):
            if isinstance(value, dict) and len(value) == 0:
                continue
            #if not isinstance(value, (str, int, float)):
            if i > 5: # 72
                del state[key]
            else:
                print(key, type(value), value)
                break
            i += 1

        #i = 0
        #for key, value in sorted(state.items()):
            #if isinstance(value, dict) and len(value) == 0:
                #continue
            #if not isinstance(value, (str, int, float)):
            #if i > 200: # 72
                #del state[key]
            #else:
                #print(key, type(value), value)
                #break
            #i += 1
        return state
