"""
Defines:
 - read_op2_geom(op2_filename=None, combine=True, subcases=None,
                 exclude_results=None, include_results=None,
                 validate=True, xref=True,
                 build_dataframe=False, skip_undefined_matrices=True,
                 mode='msc', log=None, debug=True, debug_file=None, encoding=None)
 - OP2Geom(make_geom=True, debug=False, log=None, debug_file=None, mode='msc')
   - OP2
"""
from six.moves.cPickle import dump
from pyNastran.op2.tables.geom.geom1 import GEOM1
from pyNastran.op2.tables.geom.geom2 import GEOM2
from pyNastran.op2.tables.geom.geom3 import GEOM3
from pyNastran.op2.tables.geom.geom4 import GEOM4

from pyNastran.op2.tables.geom.ept import EPT
from pyNastran.op2.tables.geom.mpt import MPT
from pyNastran.op2.tables.geom.edt import EDT
from pyNastran.op2.tables.geom.edom import EDOM

from pyNastran.op2.tables.geom.dit import DIT
from pyNastran.op2.tables.geom.dynamics import DYNAMICS

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.errors import DuplicateIDsError
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


class OP2GeomCommon(OP2, GEOM1, GEOM2, GEOM3, GEOM4, EPT, MPT, EDT, EDOM, DIT, DYNAMICS):
    """interface for the OP2Geom class for to loading subclasses"""
    def __init__(self, make_geom=True,
                 debug=False, log=None, debug_file=None, mode=None):
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
        mode : str; default=None -> 'msc'
            {msc, nx}

        """
        GEOM1.__init__(self)
        GEOM2.__init__(self)
        GEOM3.__init__(self)
        GEOM4.__init__(self)

        EPT.__init__(self)
        MPT.__init__(self)
        EDT.__init__(self)
        EDOM.__init__(self)
        DIT.__init__(self)
        DYNAMICS.__init__(self)

        OP2.__init__(self, debug=debug, log=log, debug_file=debug_file, mode=mode)
        self.make_geom = True

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


    def _read_contact_4(self, data, ndata):
        """
        reads the CONTACT/CONTACTS table
        Table of Bulk Data entry related to surface contact

        """
        return self._read_geom_4(self._contact_map, data, ndata)

    def _read_viewtb_4(self, data, ndata):
        """
        View information table
        Contains the relationship between each p-element and its view-elements and view-grids.
        """
        return self._read_geom_4(self._viewtb_map, data, ndata)


class OP2Geom(BDF, OP2GeomCommon):
    """creates an interface for the OP2 and BDF classes"""
    _properties = ['is_bdf_vectorized', 'nid_map', 'wtmass',
                   'is_real', 'is_complex', 'is_random',
                   '_sort_method', 'is_sort1', 'is_sort2',
                   'matrix_tables', 'table_name_str']
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

    @property
    def is_geometry(self):
        return True

    def read_op2(self, op2_filename=None, combine=True,
                 build_dataframe=None, skip_undefined_matrices=False, encoding=None):
        """see ``OP2.read_op2``"""
        OP2.read_op2(self, op2_filename=op2_filename, combine=combine,
                     build_dataframe=build_dataframe,
                     skip_undefined_matrices=skip_undefined_matrices,
                     encoding=encoding)
        if len(self.nodes) == 0:
            self.gpdt_to_nodes()

    def gpdt_to_nodes(self):
        """converts the GPDT & EQEXIN tables to node ids"""
        eqexin = self.op2_results.eqexin
        gpdt = self.op2_results.gpdt
        msg = ''
        if eqexin is None:
            msg += 'eqexin is None; '
        if gpdt is None:
            msg += 'gpdt is None'
            return
        if msg:
            self.log.error('Cannot convert EQEXIN/GPDT to nodes because %s' % msg.rstrip('; '))
            return

        nid_cp_cd_ps = gpdt.nid_cp_cd_ps
        xyz = gpdt.xyz
        nids = eqexin.nid
        for nid, nid_cp_cd_psi, xyzi in zip(nids, nid_cp_cd_ps, xyz):
            _nid, cp, cd, ps = nid_cp_cd_psi
            self.add_grid(nid, xyzi, cp=cp, cd=cd, ps=ps, seid=0, comment='')

    def __getstate__(self):
        """clears out a few variables in order to pickle the object"""
        raise NotImplementedError()
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        #adfasd
        #state = BDF.__getstate__(self)
        #print(state)
        #state = self.__dict__.copy()

        # Remove the unpicklable entries.
        #i = 0
        #for key, value in sorted(state.items()):
            #if isinstance(value, dict) and len(value) == 0:
                #continue
            ##if not isinstance(value, (str, int, float)):
            #if i > 5: # 72
                #del state[key]
            #else:
                #print(key, type(value), value)
                #break
            #i += 1

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
        #return state

    def export_hdf5_file(self, hdf5_file, exporter=None):
        """
        Converts the OP2 objects into hdf5 object

        Parameters
        ----------
        hdf5_file : H5File()
            an h5py object
        exporter : HDF5Exporter; default=None
            unused

        TODO: doesn't support:
          - BucklingEigenvalues

        """
        #from pyNastran.op2.op2_interface.hdf5_interface import export_op2_to_hdf5_file
        #op2_model = self
        OP2GeomCommon.export_hdf5_file(self, hdf5_file)
        BDF.export_hdf5_file(self, hdf5_file)
