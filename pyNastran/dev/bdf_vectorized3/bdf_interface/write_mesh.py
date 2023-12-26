from __future__ import annotations
from typing import Union, Optional, TYPE_CHECKING
from io import StringIO

from pyNastran.dev.bdf_vectorized3.bdf_interface.bdf_attributes import BDFAttributes
from pyNastran.bdf.bdf_interface.write_mesh import _output_helper, _fix_sizes
from pyNastran.dev.bdf_vectorized3.types import TextIOLike
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.types import TextIOLike
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


class WriteMesh(BDFAttributes):
    def __init__(self):
        super().__init__()
        #BDFAttributes.__init__(self)
        self.writer = Writer(self)

    def write_bdf(self, out_filename: Optional[Union[str, StringIO]]=None,
                  encoding: Optional[str]=None,
                  size: int=8,
                  nodes_size: Optional[int]=None,
                  elements_size: Optional[int]=None,
                  loads_size: Optional[int]=None,
                  is_double: bool=False,
                  interspersed: bool=False, enddata: Optional[bool]=None,
                  write_header: bool=True, close: bool=True) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        out_filename : varies; default=None
            str        - the name to call the output bdf
            file       - a file object
            StringIO() - a StringIO object
            None       - pops a dialog
        encoding : str; default=None -> system specified encoding
            the unicode encoding
            latin1, cp1252, and utf8 are generally good options
        size : int; {8, 16}
            the field size
        is_double : bool; default=False
            False : small field
            True : large field
        interspersed : bool; default=True
            Writes a bdf with properties & elements
            interspersed like how Patran writes the bdf.  This takes
            slightly longer than if interspersed=False, but makes it
            much easier to compare to a Patran-formatted bdf and is
            more clear.
        enddata : bool; default=None
            bool - enable/disable writing ENDDATA
            None - depends on input BDF
        write_header : bool; default=True
            flag for writing the pyNastran header
        close : bool; default=True
            should the output file be closed

        """
        is_long_ids = False
        #else:
            #is_long_ids, size = self._get_long_ids(size)

        log = self.log
        out_filename, size = _output_helper(out_filename,
                                            interspersed, size, is_double, log)
        encoding = self.get_encoding(encoding)
        #assert encoding.lower() in ['ascii', 'latin1', 'utf8'], encoding

        has_read_write = hasattr(out_filename, 'read') and hasattr(out_filename, 'write')
        if has_read_write:
            bdf_file = out_filename
        else:
            log.debug(f'---starting BDF.write_bdf of {out_filename}---')
            assert isinstance(encoding, str), encoding
            bdf_file = open(out_filename, 'w', encoding=encoding)
        self.writer._write_header(bdf_file, encoding, write_header=write_header)

        #log.warning('write_bdf')
        #-------------------------------------------------------------
        self.write_bulk_data(
            bdf_file, size=size, is_double=is_double, interspersed=interspersed,
            enddata=enddata, close=close,
            nodes_size=nodes_size, elements_size=elements_size,
            loads_size=loads_size, is_long_ids=is_long_ids)

        #if self.suport or self.suport1:
            #bdf_file.write('$CONSTRAINTS\n')
            #for suport in self.suport:
                #bdf_file.write(suport.write_card(size, is_double))
            #for unused_suport_id, suport in sorted(self.suport1.items()):
                #bdf_file.write(suport.write_card(size, is_double))
        #self._write_common(bdf_file, size=size, is_double=False, is_long_ids=None)

    def write_bulk_data(self, bdf_file: TextIOLike,
                        size: int=8, is_double: bool=False,
                        interspersed: bool=False,
                        enddata: Optional[bool]=None, close: bool=True,
                        nodes_size: Optional[int]=None,
                        elements_size: Optional[int]=None,
                        loads_size: Optional[int]=None,
                        is_long_ids: bool=False) -> None:
        """
        Writes the BDF.

        Parameters
        ----------
        bdf_file : varies
            file       - a file object
            StringIO() - a StringIO object
        size : int; {8, 16}
            the field size
        is_double : bool; default=False
            False : small field
            True : large field
        interspersed : bool; default=True
            Writes a bdf with properties & elements
            interspersed like how Patran writes the bdf.  This takes
            slightly longer than if interspersed=False, but makes it
            much easier to compare to a Patran-formatted bdf and is
            more clear.
        enddata : bool; default=None
            bool - enable/disable writing ENDDATA
            None - depends on input BDF
        close : bool; default=True
            should the output file be closed

        .. note:: is_long_ids is only needed if you have ids longer
                  than 8 characters. It's an internal parameter, but if
                  you're calling the new sub-function, you might need
                  it.  Chances are you won't.
        """
        self.writer.write_bulk_data(
            bdf_file,
            size=size, is_double=is_double,
            interspersed=interspersed,
            enddata=enddata, close=close,
            nodes_size=nodes_size,
            elements_size=elements_size,
            loads_size=loads_size,
            is_long_ids=is_long_ids)


class Writer:
    def __init__(self, model: BDF):
        self.model = model

    def write_bulk_data(self, bdf_file,
                        size: int=8, is_double: bool=False,
                        interspersed: bool=False,
                        enddata: Optional[bool]=None, close: bool=True,
                        nodes_size: Optional[int]=None,
                        elements_size: Optional[int]=None,
                        loads_size: Optional[int]=None,
                        is_long_ids: bool=False) -> None:
        size, nodes_size, elements_size, loads_size = _fix_sizes(
            size, nodes_size, elements_size, loads_size)

        model = self.model

        if model.run_testing_checks:
            model.quality()
            #model.spcadd.get_reduced_spcs()

        self._write_params(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_nodes(bdf_file, nodes_size, is_double, is_long_ids=is_long_ids)
        self._write_elements(bdf_file, elements_size, is_double, is_long_ids=is_long_ids)
        self._write_properties(bdf_file, size, is_double, is_long_ids=is_long_ids)
        #self._write_properties_by_element_type(bdf_file, size, is_double, is_long_ids)
        self._write_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_masses(bdf_file, size, is_double, is_long_ids=is_long_ids)

        # split out for write_bdf_symmetric
        self._write_rigid_elements(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_nonstructural_mass(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_aero(bdf_file, size, is_double, is_long_ids=is_long_ids)

        self._write_common(bdf_file, loads_size, is_double, is_long_ids=is_long_ids)
        if (enddata is None and 'ENDDATA' in model.card_count) or enddata:
            bdf_file.write('ENDDATA\n')
        if close:
            bdf_file.close()

    def _write_header(self, bdf_file: TextIOLike,
                      encoding: str, write_header: bool=True) -> None:
        """Writes the executive and case control decks."""
        self._set_punch()

        model = self.model
        if model.nastran_format and write_header:
            bdf_file.write(f'$pyNastran: version={model.nastran_format}\n')
            bdf_file.write(f'$pyNastran: punch={model.punch}\n')
            bdf_file.write(f'$pyNastran: encoding={encoding}\n')
            #bdf_file.write(f'$pyNastran: nnodes={model.grid.n:d}\n')
            #bdf_file.write(f'$pyNastran: nelements={len(self.elements):d}\n')

        if not model.punch:
            self._write_executive_control_deck(bdf_file)
            self._write_case_control_deck(bdf_file)

    def _write_common(self, bdf_file: TextIOLike,
                      size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """
        Write the common outputs so none get missed...

        Parameters
        ----------
        bdf_file : file
            the file object
        size : int (default=8)
            the field width
        is_double : bool (default=False)
            is this double precision

        """
        #model = self.model
        self._write_dmigs(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_loads(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_dynamic(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_aero_control(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_static_aero(bdf_file, size, is_double, is_long_ids=is_long_ids)

        write_aero_in_flutter, write_aero_in_gust = self.find_aero_location()
        self._write_flutter(bdf_file, size, is_double, write_aero_in_flutter,
                            is_long_ids=is_long_ids)
        self._write_gust(bdf_file, size, is_double, write_aero_in_gust, is_long_ids=is_long_ids)

        self._write_thermal(bdf_file, size, is_double, is_long_ids=is_long_ids)
        #self._write_thermal_materials(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_constraints(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_optimization(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_tables(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_sets(bdf_file, size, is_double, is_long_ids=is_long_ids)
        #self._write_superelements(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_contact(bdf_file, size, is_double, is_long_ids=is_long_ids)
        #self._write_parametric(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_rejects(bdf_file, size, is_double, is_long_ids=is_long_ids)
        self._write_coords(bdf_file, size, is_double, is_long_ids=is_long_ids)

        #if self.acmodl:
            #bdf_file.write(self.acmodl.write_card(size, is_double))

    def find_aero_location(self) -> tuple[bool, bool]:
        """Determines where the AERO card should be written"""
        model = self.model
        write_aero_in_flutter = False
        write_aero_in_gust = False
        if model.aero:
            if len(model.flfact) or len(model.flutter) or model.mkaeros:
                write_aero_in_flutter = True
            elif len(model.gust):
                write_aero_in_gust = True
            else:
                # an AERO card exists, but no FLUTTER, FLFACT, MKAEROx or GUST card
                write_aero_in_flutter = True
        return write_aero_in_flutter, write_aero_in_gust

    def _write_executive_control_deck(self, bdf_file: TextIOLike) -> None:
        """Writes the executive control deck."""
        msg = ''
        model = self.model
        for line in model.system_command_lines:
            msg += line + '\n'

        if model.executive_control_lines:
            msg += '$EXECUTIVE CONTROL DECK\n'

            if model.sol_iline is not None:
                if model.sol == 600 and model.sol_method:
                    new_sol = f'SOL 600,{model.sol_method}'
                else:
                    new_sol = f'SOL {model.sol}'
                model.executive_control_lines[model.sol_iline] = new_sol

            for line in model.executive_control_lines:
                msg += line + '\n'
        elif model.sol is not None:
            msg += '$EXECUTIVE CONTROL DECK\n'
            msg += f'SOL {model.sol}\n'
            msg += 'CEND\n'
        bdf_file.write(msg)

    def _write_case_control_deck(self, bdf_file: TextIOLike) -> None:
        """Writes the Case Control Deck."""
        model = self.model
        case_control_deck = model.case_control_deck
        if case_control_deck:
            msg = '$CASE CONTROL DECK\n'
            if model.superelement_models:
                msg += case_control_deck.write(write_begin_bulk=False)
            else:
                msg += str(case_control_deck)
                assert 'BEGIN BULK' in msg, msg
            bdf_file.write(''.join(msg))

    def _set_punch(self) -> None:
        """updates the punch flag"""
        model = self.model
        if model.punch is None:
            # writing a mesh without using read_bdf
            if model.system_command_lines or model.executive_control_lines or model.case_control_deck:
                model.punch = False
            else:
                model.punch = True

    def _write_params(self, bdf_file: TextIOLike,
                      size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the PARAM cards"""
        model = self.model
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        if model.params or model.dti or model.mdlprm:
            bdf_file.write('$PARAMS\n')
            for unused_name, dti in sorted(model.dti.items()):
                bdf_file.write(dti.write_card(size=size, is_double=is_double))

            for (unused_key, param) in sorted(model.params.items()):
                bdf_file.write(param.write_card(size, is_double))
            if model.mdlprm:
                bdf_file.write(model.mdlprm.write_card(size, is_double))

    def _write_nodes(self, bdf_file: TextIOLike,
                     size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the NODE-type cards"""
        model = self.model
        if model.spoint.n:
            bdf_file.write('$SPOINTS\n')
            model.spoint.write_file(bdf_file, size=size)
        if model.epoint.n:
            bdf_file.write('$EPOINTS\n')
            model.epoint.write_file(bdf_file, size=size)
        if model.point.n:
            bdf_file.write('$POINTS\n')
            model.point.write_file(bdf_file, size=size)

        #if self._is_axis_symmetric:
            #if self.axic:
                #bdf_file.write(self.axic.write_card(size, is_double))
            #if self.axif:
                #bdf_file.write(self.axif.write_card(size, is_double))
            #for unused_nid, ringax_pointax in sorted(self.ringaxs.items()):
                #bdf_file.write(ringax_pointax.write_card(size, is_double))
            #for unused_ringfl, ringfl in sorted(self.ringfl.items()):
                #bdf_file.write(ringfl.write_card(size, is_double))
            #for unused_nid, gridb in sorted(self.gridb.items()):
                #bdf_file.write(gridb.write_card(size, is_double))
        #if self.cyax:
            #bdf_file.write(self.cyax.write_card(size, is_double))

        self._write_grids(bdf_file, size=size, is_double=is_double)
        #if self.seqgp:
            #bdf_file.write(self.seqgp.write_card(size, is_double))

    def _write_grids(self, bdf_file: TextIOLike,
                     size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the GRID-type cards"""
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if model.grid.n:
            bdf_file.write('$NODES\n')
            #if self.grdset:
                #bdf_file.write(self.grdset.write_card(size))
            model.grid.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

    def _write_elements(self, bdf_file: TextIOLike,
                        size: int=8, is_double: bool=False,
                        is_long_ids: Optional[bool]=None) -> None:
        """Writes the elements in a sorted order"""
        model = self.model
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        # plot
        model.plotel.write_file(bdf_file, size=size, is_double=is_double)
        model.plotel3.write_file(bdf_file, size=size, is_double=is_double)
        model.plotel4.write_file(bdf_file, size=size, is_double=is_double)
        model.plotel6.write_file(bdf_file, size=size, is_double=is_double)
        model.plotel8.write_file(bdf_file, size=size, is_double=is_double)

        # celas
        model.celas1.write_file(bdf_file, size=size, is_double=is_double)
        model.celas2.write_file(bdf_file, size=size, is_double=is_double)
        model.celas3.write_file(bdf_file, size=size, is_double=is_double)
        model.celas4.write_file(bdf_file, size=size, is_double=is_double)

        # cdamp
        model.cdamp1.write_file(bdf_file, size=size, is_double=is_double)
        model.cdamp2.write_file(bdf_file, size=size, is_double=is_double)
        model.cdamp3.write_file(bdf_file, size=size, is_double=is_double)
        model.cdamp4.write_file(bdf_file, size=size, is_double=is_double)
        model.cdamp5.write_file(bdf_file, size=size, is_double=is_double)

        model.cvisc.write_file(bdf_file, size=size, is_double=is_double)
        model.cgap.write_file(bdf_file, size=size, is_double=is_double)

        # bush
        model.cbush.write_file(bdf_file, size=size, is_double=is_double)
        model.cbush1d.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.cbush2d.write(size=size))

        # fast
        model.cfast.write_file(bdf_file, size=size, is_double=is_double)

        # rod
        model.crod.write_file(bdf_file, size=size, is_double=is_double)
        model.conrod.write_file(bdf_file, size=size, is_double=is_double)
        model.ctube.write_file(bdf_file, size=size, is_double=is_double)

        # bar/beam/shear
        model.cbar.write_file(bdf_file, size=size, is_double=is_double)
        model.cbarao.write_file(bdf_file, size=size, is_double=is_double)
        model.cbeam.write_file(bdf_file, size=size, is_double=is_double)
        model.cbend.write_file(bdf_file, size=size, is_double=is_double)
        model.cshear.write_file(bdf_file, size=size, is_double=is_double)

        # shells
        #all_shells = (mo
        model.ctria3.write_file(bdf_file, size=size, is_double=is_double)
        model.cquad4.write_file(bdf_file, size=size, is_double=is_double)
        model.ctria6.write_file(bdf_file, size=size, is_double=is_double)
        model.cquad8.write_file(bdf_file, size=size, is_double=is_double)
        model.ctriar.write_file(bdf_file, size=size, is_double=is_double)
        model.cquadr.write_file(bdf_file, size=size, is_double=is_double)
        model.cquad.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.snorm.write(size=size), is_double=is_double)

        # axisymmetric shells
        model.ctriax.write_file(bdf_file, size=size, is_double=is_double)
        model.ctriax6.write_file(bdf_file, size=size, is_double=is_double)
        model.cquadx.write_file(bdf_file, size=size, is_double=is_double)
        model.cquadx4.write_file(bdf_file, size=size, is_double=is_double)
        model.cquadx8.write_file(bdf_file, size=size, is_double=is_double)
        model.ctrax3.write_file(bdf_file, size=size, is_double=is_double)
        model.ctrax6.write_file(bdf_file, size=size, is_double=is_double)

        # plate stress
        model.cplsts3.write_file(bdf_file, size=size, is_double=is_double)
        model.cplsts4.write_file(bdf_file, size=size, is_double=is_double)
        model.cplsts6.write_file(bdf_file, size=size, is_double=is_double)
        model.cplsts8.write_file(bdf_file, size=size, is_double=is_double)

        # plate strain
        model.cplstn3.write_file(bdf_file, size=size, is_double=is_double)
        model.cplstn4.write_file(bdf_file, size=size, is_double=is_double)
        model.cplstn6.write_file(bdf_file, size=size, is_double=is_double)
        model.cplstn8.write_file(bdf_file, size=size, is_double=is_double)

        # solids
        model.ctetra.write_file(bdf_file, size=size, is_double=is_double)
        model.cpenta.write_file(bdf_file, size=size, is_double=is_double)
        model.chexa.write_file(bdf_file, size=size, is_double=is_double)
        model.cpyram.write_file(bdf_file, size=size, is_double=is_double)

        # nx solid cohesive zone
        model.cpentcz.write_file(bdf_file, size=size, is_double=is_double)
        model.chexcz.write_file(bdf_file, size=size, is_double=is_double)
        # msc solid cohesive zone
        model.cifpent.write_file(bdf_file, size=size, is_double=is_double)
        model.cifhex.write_file(bdf_file, size=size, is_double=is_double)

        # acoustic shells
        model.caabsf.write_file(bdf_file, size=size, is_double=is_double)
        model.paabsf.write_file(bdf_file, size=size, is_double=is_double)

        # acoustic solids
        model.chacab.write_file(bdf_file, size=size, is_double=is_double)
        model.chacbr.write_file(bdf_file, size=size, is_double=is_double)

        # other
        model.genel.write_file(bdf_file, size=size, is_double=is_double)

    def _write_rejects(self, bdf_file: TextIOLike,
                       size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the rejected (processed) cards and the rejected unprocessed
        cardlines

        """
        from pyNastran.bdf.field_writer import print_card_8, print_card_16
        if size == 8:
            print_func = print_card_8
        else:
            print_func = print_card_16

        model = self.model
        if model.reject_cards:
            bdf_file.write('$REJECT_CARDS\n')
            for reject_card in model.reject_cards:
                try:
                    bdf_file.write(print_func(reject_card))
                except RuntimeError:
                    if len(reject_card) > 0:
                        line0 = reject_card[0].upper()
                        if line0.startswith('ADAPT'):
                            for line in reject_card:
                                assert isinstance(line, str), line
                                bdf_file.write(line+'\n')
                            continue
                    for field in reject_card:
                        if field is not None and '=' in field:
                            raise SyntaxError('cannot reject equal signed '
                                              'cards\ncard=%s\n' % reject_card)
                    raise

        if model.reject_lines:
            bdf_file.write('$REJECT_LINES\n')
            for reject_lines in model.reject_lines:
                if isinstance(reject_lines, (list, tuple)):
                    for reject in reject_lines:
                        reject2 = reject.rstrip()
                        if reject2:
                            bdf_file.write('%s\n' % reject2)
                elif isinstance(reject_lines, str):
                    reject2 = reject_lines.rstrip()
                    if reject2:
                        bdf_file.write('%s\n' % reject2)
                else:
                    raise TypeError(reject_lines)

    def _write_rigid_elements(self, bdf_file: TextIOLike,
                              size: int=8, is_double: bool=False,
                              is_long_ids: Optional[bool]=None) -> None:
        """Writes the rigid elements in a sorted order"""
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        model.rrod.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.rrod1.write(size=size))
        model.rbar.write_file(bdf_file, size=size, is_double=is_double)
        model.rbar1.write_file(bdf_file, size=size, is_double=is_double)
        model.rbe1.write_file(bdf_file, size=size, is_double=is_double)
        model.rbe2.write_file(bdf_file, size=size, is_double=is_double)
        model.rbe3.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.rspline.write(size=size))
        #bdf_file.write(model.rsscon.write(size=size))

    def _write_nonstructural_mass(self, bdf_file: TextIOLike,
                                  size: int=8, is_double: bool=False,
                                  is_long_ids: Optional[bool]=None) -> None:
        """Writes the rigid elements in a sorted order"""
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        model.nsmadd.write_file(bdf_file, size=size, is_double=is_double)
        model.nsm.write_file(bdf_file, size=size, is_double=is_double)
        model.nsm1.write_file(bdf_file, size=size, is_double=is_double)
        model.nsml.write_file(bdf_file, size=size, is_double=is_double)
        model.nsml1.write_file(bdf_file, size=size, is_double=is_double)

    def _write_masses(self, bdf_file: TextIOLike,
                      size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the mass cards sorted by ID"""
        model = self.model
        all_masses = (
            model.conm1, model.conm2,
            model.pmass,
            model.cmass1, model.cmass2, model.cmass3, model.cmass4,
        )
        masses = [mass for mass in all_masses if mass.n > 0]
        if len(masses) == 0:
            return
        bdf_file.write('$MASSES\n')
        for mass in masses:
            mass.write_file(bdf_file, size=size, is_double=is_double)

        #bdf_file.write(model.conm1.write(size=size))
        #bdf_file.write(model.conm2.write(size=size))
        #bdf_file.write(model.pmass.write(size=size))
        #bdf_file.write(model.cmass1.write(size=size))
        #bdf_file.write(model.cmass2.write(size=size))
        #bdf_file.write(model.cmass3.write(size=size))
        #bdf_file.write(model.cmass4.write(size=size))
        return
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

    def _write_properties(self, bdf_file: TextIOLike,
                          size: int=8, is_double: bool=False,
                          is_long_ids: Optional[bool]=None) -> None:
        """Writes the properties in a sorted order"""
        model = self.model
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        # spring/damp
        model.pelas.write_file(bdf_file, size=size, is_double=is_double)
        model.pelast.write_file(bdf_file, size=size, is_double=is_double)
        model.pdamp.write_file(bdf_file, size=size, is_double=is_double)
        model.pdampt.write_file(bdf_file, size=size, is_double=is_double)
        model.pvisc.write_file(bdf_file, size=size, is_double=is_double)
        model.pgap.write_file(bdf_file, size=size, is_double=is_double)

        # bush
        model.pbush.write_file(bdf_file, size=size, is_double=is_double)
        model.pbusht.write_file(bdf_file, size=size, is_double=is_double)
        model.pbush1d.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.pbush2d.write(size=size))

        # fast
        model.pfast.write_file(bdf_file, size=size, is_double=is_double)

        # rod
        model.prod.write_file(bdf_file, size=size, is_double=is_double)
        model.ptube.write_file(bdf_file, size=size, is_double=is_double)

        # bar
        model.pbar.write_file(bdf_file, size=size, is_double=is_double)
        model.pbarl.write_file(bdf_file, size=size, is_double=is_double)
        model.pbrsect.write_file(bdf_file, size=size, is_double=is_double)

        # beam
        model.pbeam.write_file(bdf_file, size=size, is_double=is_double)
        model.pbeaml.write_file(bdf_file, size=size, is_double=is_double)
        model.pbcomp.write_file(bdf_file, size=size, is_double=is_double)

        # bend
        model.pbend.write_file(bdf_file, size=size, is_double=is_double)

        # shear
        model.pshear.write_file(bdf_file, size=size, is_double=is_double)

        # shell
        model.pshell.write_file(bdf_file, size=size, is_double=is_double)
        model.pcomp.write_file(bdf_file, size=size, is_double=is_double)
        model.pcompg.write_file(bdf_file, size=size, is_double=is_double)

        # planar shells
        model.plplane.write_file(bdf_file, size=size, is_double=is_double)
        model.pplane.write_file(bdf_file, size=size, is_double=is_double)
        model.pshln1.write_file(bdf_file, size=size, is_double=is_double)
        model.pshln2.write_file(bdf_file, size=size, is_double=is_double)

        # solid
        model.psolid.write_file(bdf_file, size=size, is_double=is_double)
        model.plsolid.write_file(bdf_file, size=size, is_double=is_double)
        model.psolcz.write_file(bdf_file, size=size, is_double=is_double)
        model.pcomps.write_file(bdf_file, size=size, is_double=is_double)
        model.pcompls.write_file(bdf_file, size=size, is_double=is_double)

        # acoustic
        model.pacbar.write_file(bdf_file, size=size, is_double=is_double)

    def _write_materials(self, bdf_file: TextIOLike, size: int=8, is_double: bool=False,
                         is_long_ids: Optional[bool]=None) -> None:
        """Writes the materials in a sorted order"""
        model = self.model
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)

        structural_materials = [
            model.mat1, model.mat2, model.mat3,
            model.mat8, model.mat9, model.mat10, model.mat11,
            model.mat10c, model.matort,
            model.mathe, # hyperelastic
            model.mathp,  # hyperelastic
        ]
        thermal_materials = [model.mat4, model.mat5] + model.thermal_material_dep_cards
        is_materials = any([mat.n for mat in structural_materials])
        is_thermal_materials = any([mat.n for mat in thermal_materials])
        if is_materials:
            bdf_file.write('$MATERIALS\n')
            model.mat1.write_file(bdf_file, size=size, is_double=is_double)  # isotropic
            model.mat2.write_file(bdf_file, size=size, is_double=is_double)
            model.mat3.write_file(bdf_file, size=size, is_double=is_double)
            model.mat8.write_file(bdf_file, size=size, is_double=is_double)  # orthotropic
            model.mat9.write_file(bdf_file, size=size, is_double=is_double)
            model.mat10.write_file(bdf_file, size=size, is_double=is_double)
            model.mat11.write_file(bdf_file, size=size, is_double=is_double)
            model.mat10c.write_file(bdf_file, size=size, is_double=is_double)
            model.matort.write_file(bdf_file, size=size, is_double=is_double)
             #hyperelastic
            model.mathe.write_file(bdf_file, size=size, is_double=is_double)
            model.mathp.write_file(bdf_file, size=size, is_double=is_double)

        if is_thermal_materials:
            bdf_file.write('$THERMAL_MATERIALS\n')
            model.mat4.write_file(bdf_file, size=size, is_double=is_double)
            model.mat5.write_file(bdf_file, size=size, is_double=is_double)

        model.mats1.write_file(bdf_file, size=size, is_double=is_double)

        model.matt1.write_file(bdf_file, size=size, is_double=is_double)
        model.matt2.write_file(bdf_file, size=size, is_double=is_double)
        model.matt3.write_file(bdf_file, size=size, is_double=is_double)
        model.matt4.write_file(bdf_file, size=size, is_double=is_double)
        model.matt5.write_file(bdf_file, size=size, is_double=is_double)
        model.matt8.write_file(bdf_file, size=size, is_double=is_double)
        model.matt9.write_file(bdf_file, size=size, is_double=is_double)

        #is_big_materials = hasattr(model, 'big_materials') and model.big_materials
        is_materials_dict = len(model.nxstrats)
        #is_materials = (self.materials or self.hyperelastic_materials or self.creep_materials or
                        #self.MATS3 or self.MATS8 or
                        #self.MATT3 or self.MATT4 or self.MATT5 or
                        #self.nxstrats or is_big_materials)

        if is_materials_dict:
            #for (unused_mid, material) in sorted(self.hyperelastic_materials.items()):
                #bdf_file.write(material.write_card(size, is_double))
            #for (unused_mid, material) in sorted(self.creep_materials.items()):
                #bdf_file.write(material.write_card(size, is_double))

            #for (unused_mid, material) in sorted(self.MATS3.items()):
                #bdf_file.write(material.write_card(size, is_double))
            #for (unused_mid, material) in sorted(self.MATS8.items()):
                #bdf_file.write(material.write_card(size, is_double))

            #for (unused_mid, material) in sorted(self.MATT3.items()):
                #bdf_file.write(material.write_card(size, is_double))
            #for (unused_mid, material) in sorted(self.MATT4.items()):
                #bdf_file.write(material.write_card(size, is_double))
            #for (unused_mid, material) in sorted(self.MATT5.items()):
                #bdf_file.write(material.write_card(size, is_double))
            for (unused_sid, nxstrat) in sorted(model.nxstrats.items()):
                bdf_file.write(nxstrat.write_card(size, is_double))

            #if is_big_materials:
                #for unused_mid, mat in sorted(self.big_materials.items()):
                    #bdf_file.write(mat.write_card_16(is_double))

    def _write_sets(self, bdf_file: TextIOLike, size: int=8, is_double: bool=False,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the SETx cards sorted by ID"""
        model = self.model
        model.aset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.bset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.cset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.qset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.omit.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.uset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.set1.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        #bdf_file.write(model.set2.write(size=size))  #  faked
        model.set3.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.seset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.sebset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.secset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.seqset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        #model.seuset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.release.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

    def _write_superelements(self, bdf_file: TextIOLike, size: int=8, is_double: bool=False,
                             is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the Superelement cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        pass

    def _write_tables(self, bdf_file: TextIOLike, size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the TABLEx cards sorted by ID"""
        model = self.model
        if model.tables or model.tables_d or model.tables_m or model.tables_sdamping:
            bdf_file.write('$TABLES\n')
            for (unused_id, table) in sorted(model.tables.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(model.tables_d.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(model.tables_m.items()):
                bdf_file.write(table.write_card(size, is_double))
            for (unused_id, table) in sorted(model.tables_sdamping.items()):
                bdf_file.write(table.write_card(size, is_double))

        #if model.random_tables:
            #bdf_file.write('$RANDOM TABLES\n')
            #for (unused_id, table) in sorted(model.random_tables.items()):
                #bdf_file.write(table.write_card(size, is_double))

    def _write_thermal(self, bdf_file: TextIOLike, size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal cards"""
        model = self.model
        is_thermal = any([card.n > 0 for card in model.thermal_element_cards])
        if is_thermal:
            bdf_file.write('$THERMAL_ELEMENTS\n')
            model.chbdye.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.chbdyp.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.chbdyg.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.phbdy.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        is_thermal = any([card.n > 0 for card in model.thermal_boundary_condition_cards])
        if is_thermal:
            model.conv.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.pconv.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.convm.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.pconvm.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.radcav.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.radlst.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        #model.radmtx.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.radbc.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.radm.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.radset.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.tempbc.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

        model.view.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
        model.view3d.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

    def _write_thermal_materials(self, bdf_file: TextIOLike,
                                 size: int=8, is_double: bool=False,
                                 is_long_ids: Optional[bool]=None) -> None:
        """Writes the thermal materials in a sorted order"""
        pass

    def _write_aero(self, bdf_file: TextIOLike,
                    size: int=8, is_double: bool=False,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero cards"""
        model = self.model
        aero_cards = model.monitor_point_cards

        is_aero = any([card.n > 0 for card in aero_cards])
        if is_aero:
            bdf_file.write('$AERO\n')
            model.monpnt1.write_file(bdf_file, size=size)
            model.monpnt2.write_file(bdf_file, size=size)
            model.monpnt3.write_file(bdf_file, size=size)

        #if model.paeros: # or model.caeros or model.splines:
            #for (unused_id, caero) in sorted(model.caeros.items()):
                #bdf_file.write(caero.write_card(size, is_double))
            #for (unused_id, spline) in sorted(model.splines.items()):
                #bdf_file.write(spline.write_card(size, is_double))
        return
        self.zona.write_bdf(bdf_file, size=8, is_double=False)

    def _write_aero_control(self, bdf_file: TextIOLike,
                            size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the aero control surface cards"""
        model = self.model
        if(len(model.aefact) or len(model.aeparm) or len(model.aelink) or
           len(model.aestat) or len(model.aesurf) or len(model.aesurfs)):
            bdf_file.write('$AERO CONTROL SURFACES\n')

            #for (unused_id, aecomp) in sorted(model.aecomps.items()):
                #bdf_file.write(aecomp.write_card(size, is_double))

            model.aeparm.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.aestat.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.aesurf.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.aesurfs.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.aefact.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)
            model.aelink.write_file(bdf_file, size=size, is_double=is_double, write_card_header=False)

    def _write_static_aero(self, bdf_file: TextIOLike,
                           size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the static aero cards"""
        model = self.model
        if model.aeros or len(model.trim) or len(model.diverg) or len(model.csschd):
            bdf_file.write('$STATIC AERO\n')
            model.trim.write_file(bdf_file, size=size, is_double=is_double)
            #model.trim2.write_file(bdf_file, size=size, is_double=is_double)
            model.csschd.write_file(bdf_file, size=size, is_double=is_double)
            model.diverg.write_file(bdf_file, size=size, is_double=is_double)

            # static aero
            if model.aeros:
                bdf_file.write(model.aeros.write_card(size, is_double))

    def _write_flutter(self, bdf_file: TextIOLike,
                       size: int=8, is_double: bool=False,
                       write_aero_in_flutter: bool=True,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the flutter cards"""
        model = self.model

        model.paero1.write_file(bdf_file, size=size, is_double=is_double)
        model.paero2.write_file(bdf_file, size=size, is_double=is_double)
        model.paero3.write_file(bdf_file, size=size, is_double=is_double)
        model.paero4.write_file(bdf_file, size=size, is_double=is_double)
        model.paero5.write_file(bdf_file, size=size, is_double=is_double)

        model.caero1.write_file(bdf_file, size=size, is_double=is_double)
        model.caero2.write_file(bdf_file, size=size, is_double=is_double)
        model.caero3.write_file(bdf_file, size=size, is_double=is_double)
        model.caero4.write_file(bdf_file, size=size, is_double=is_double)
        model.caero5.write_file(bdf_file, size=size, is_double=is_double)
        model.caero7.write_file(bdf_file, size=size, is_double=is_double)  # zona

        model.spline1.write_file(bdf_file, size=size, is_double=is_double)
        model.spline2.write_file(bdf_file, size=size, is_double=is_double)
        model.spline3.write_file(bdf_file, size=size, is_double=is_double)
        model.spline4.write_file(bdf_file, size=size, is_double=is_double)
        model.spline5.write_file(bdf_file, size=size, is_double=is_double)
        #bdf_file.write(model.spline6.write(size=size))

        model.aecomp.write_file(bdf_file, size=size, is_double=is_double)  # helper for AECOMP
        model.aecompl.write_file(bdf_file, size=size, is_double=is_double)  # helper for AECOMPL
        model.aelist.write_file(bdf_file, size=size, is_double=is_double)  # aero boxes for AESURF
        model.flfact.write_file(bdf_file, size=size, is_double=is_double)  # Mach, vel, rho for FLUTTER
        model.aefact.write_file(bdf_file, size=size, is_double=is_double)

        if (write_aero_in_flutter and model.aero) or len(model.flutter) or model.mkaeros:
            bdf_file.write('$FLUTTER\n')
            if write_aero_in_flutter:
                bdf_file.write(model.aero.write_card(size, is_double))
            model.flutter.write_file(bdf_file, size=size, is_double=is_double)
            for mkaero in model.mkaeros:
                bdf_file.write(mkaero.write_card(size, is_double))

    def _write_gust(self, bdf_file: TextIOLike,
                    size: int=8, is_double: bool=False,
                    write_aero_in_gust: bool=True,
                    is_long_ids: Optional[bool]=None) -> None:
        """Writes the gust cards"""
        model = self.model
        if (write_aero_in_gust and model.aero) or len(model.gust):
            bdf_file.write('$GUST\n')
            if write_aero_in_gust:
                if model.aero is not None:
                    bdf_file.write(model.aero.write_card(size, is_double))
            bdf_file.write(model.gust.write(size=size))

    def _write_loads(self, bdf_file: TextIOLike,
                     size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """Writes the load cards sorted by ID"""
        model = self.model
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        #if self.load_combinations or self.loads or self.cyjoin:
        if any(load.n for load in model.load_cards):
            bdf_file.write('$LOADS\n')
            model.cload.write_file(bdf_file, size=size, is_double=is_double)
            model.load.write_file(bdf_file, size=size, is_double=is_double)
            model.grav.write_file(bdf_file, size=size, is_double=is_double)
            model.accel.write_file(bdf_file, size=size, is_double=is_double)
            model.accel1.write_file(bdf_file, size=size, is_double=is_double)

            model.force.write_file(bdf_file, size=size, is_double=is_double)
            model.force1.write_file(bdf_file, size=size, is_double=is_double)
            model.force2.write_file(bdf_file, size=size, is_double=is_double)
            model.moment.write_file(bdf_file, size=size, is_double=is_double)
            model.moment1.write_file(bdf_file, size=size, is_double=is_double)
            model.moment2.write_file(bdf_file, size=size, is_double=is_double)
            model.pload.write_file(bdf_file, size=size, is_double=is_double)
            model.pload1.write_file(bdf_file, size=size, is_double=is_double)
            model.pload2.write_file(bdf_file, size=size, is_double=is_double)
            model.pload4.write_file(bdf_file, size=size, is_double=is_double)
            model.sload.write_file(bdf_file, size=size, is_double=is_double)
            model.spcd.write_file(bdf_file, size=size, is_double=is_double)
            model.deform.write_file(bdf_file, size=size, is_double=is_double)

            # axisymmetric loads
            #bdf_file.write(model.ploadx1.write(size=size))

            # static thermal loads
            model.tempd.write_file(bdf_file, size=size, is_double=is_double)  # default temp
            model.temp.write_file(bdf_file, size=size, is_double=is_double)

            # dynamic thermal loads
            #bdf_file.write(model.dtemp.write(size=size))  # has nodes
            model.qhbdy.write_file(bdf_file, size=size, is_double=is_double)
            model.qbdy1.write_file(bdf_file, size=size, is_double=is_double)
            model.qbdy2.write_file(bdf_file, size=size, is_double=is_double)
            model.qbdy3.write_file(bdf_file, size=size, is_double=is_double)
            model.qvol.write_file(bdf_file, size=size, is_double=is_double)
            model.qvect.write_file(bdf_file, size=size, is_double=is_double)

            # static load sequence
            model.lseq.write_file(bdf_file, size=size, is_double=is_double)

            # rotational forces - static
            model.rforce.write_file(bdf_file, size=size, is_double=is_double)
            model.rforce1.write_file(bdf_file, size=size, is_double=is_double)
            #model.rforce2.write_file(bdf_file, size=size, is_double=is_double)

        if any(card.n for card in model.dynamic_cards):
            model.tic.write_file(bdf_file, size=size, is_double=is_double)
            model.tf.write_file(bdf_file, size=size, is_double=is_double)
            model.dphase.write_file(bdf_file, size=size, is_double=is_double)
            model.delay.write_file(bdf_file, size=size, is_double=is_double)

        if any(load.n for load in model.dynamic_load_cards):
            model.dload.write_file(bdf_file, size=size, is_double=is_double)
            model.darea.write_file(bdf_file, size=size, is_double=is_double)
            model.tload1.write_file(bdf_file, size=size, is_double=is_double)
            model.tload2.write_file(bdf_file, size=size, is_double=is_double)
            #model.tload3.write_file(bdf_file, size=size, is_double=is_double)
            model.rload1.write_file(bdf_file, size=size, is_double=is_double)
            model.rload2.write_file(bdf_file, size=size, is_double=is_double)
            #model.rloadex.write_file(bdf_file, size=size, is_double=is_double)

            # random loads
            model.randps.write_file(bdf_file, size=size, is_double=is_double)
            model.acsrce.write_file(bdf_file, size=size, is_double=is_double)

            #for (key, load_combinations) in sorted(self.load_combinations.items()):
                #for load_combination in load_combinations:
                    #try:
                        #bdf_file.write(load_combination.write_card(size, is_double))
                    #except Exception:
                        #print(f'failed printing load...type={load_combination.type} key={key!r}')
                        #raise

            #if is_long_ids:
                #for (key, loadcase) in sorted(self.loads.items()):
                    #for load in loadcase:
                        #try:
                            #bdf_file.write(load.write_card_16(is_double))
                        #except Exception:
                            #print(f'failed printing load...type={load.type} key={key!r}')
                            #raise
            #else:
                #for (key, loadcase) in sorted(self.loads.items()):
                    #for load in loadcase:
                        #try:
                            #bdf_file.write(load.write_card(size, is_double))
                        #except Exception:
                            #print(f'failed printing load...type={load.type} key={key!r}')
                            #raise

            #for unused_key, cyjoin in sorted(self.cyjoin.items()):
                #bdf_file.write(cyjoin.write_card(size, is_double))

    def _write_constraints(self, bdf_file: TextIOLike,
                           size: int=8, is_double: bool=False,
                           is_long_ids: Optional[bool]=None) -> None:
        """Writes the constraint cards sorted by ID"""
        #size, is_long_ids = self._write_mesh_long_ids_size(size, is_long_ids)
        model = self.model
        if any((card.n for card in model.spc_cards)):
            bdf_file.write('$SPCs\n')
            model.spcadd.write_file(bdf_file, size=size, is_double=is_double)
            model.spc.write_file(bdf_file, size=size, is_double=is_double)
            model.spc1.write_file(bdf_file, size=size, is_double=is_double)
            model.spcoff.write_file(bdf_file, size=size, is_double=is_double)
            model.bndfix.write_file(bdf_file, size=size, is_double=is_double)
            model.bndfree.write_file(bdf_file, size=size, is_double=is_double)

        if any((card.n for card in model.mpc_cards)):
            bdf_file.write('$MPCs\n')
            model.mpcadd.write_file(bdf_file, size=size, is_double=is_double)
            model.mpc.write_file(bdf_file, size=size, is_double=is_double)
            #bdf_file.write(model.mpcax.write(size=size))

        if len(model.suport):
            bdf_file.write('$CONSTRAINTS\n')
            model.suport.write_file(bdf_file, size=size, is_double=is_double)

    def _write_contact(self, bdf_file: TextIOLike,
                       size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the contact cards sorted by ID"""
        model = self.model
        is_contacta = max([card.n for card in model.contact_cards])
        is_contactb = model.bctparm or model.bctpara
        # or self.bcparas
        if is_contacta or is_contactb:
            bdf_file.write('$CONTACT\n')
            # msc contact
            model.bcbody.write_file(bdf_file, size=size, is_double=is_double)
            model.bcbody1.write_file(bdf_file, size=size, is_double=is_double)

            model.bcrpara.write_file(bdf_file, size=size, is_double=is_double)
            for (unused_id, bctparam) in sorted(model.bctparm.items()):
                bdf_file.write(bctparam.write_card(size, is_double))
            for (unused_id, bctpara) in sorted(model.bctpara.items()):
                bdf_file.write(bctpara.write_card(size, is_double))

            # general contact
            model.bctset.write_file(bdf_file, size=size, is_double=is_double)
            model.bctadd.write_file(bdf_file, size=size, is_double=is_double)

            # glue contact
            model.bedge.write_file(bdf_file, size=size, is_double=is_double)
            model.bgset.write_file(bdf_file, size=size, is_double=is_double)
            model.bgadd.write_file(bdf_file, size=size, is_double=is_double)

            # surfaces/solids
            model.bsurf.write_file(bdf_file, size=size, is_double=is_double)
            model.bsurfs.write_file(bdf_file, size=size, is_double=is_double)
            model.bcprop.write_file(bdf_file, size=size, is_double=is_double)
            model.bcprops.write_file(bdf_file, size=size, is_double=is_double)
            model.bconp.write_file(bdf_file, size=size, is_double=is_double)
            model.blseg.write_file(bdf_file, size=size, is_double=is_double)
            model.bfric.write_file(bdf_file, size=size, is_double=is_double)

    def _write_coords(self, bdf_file: TextIOLike,
                      size: int=8, is_double: bool=False,
                      is_long_ids: Optional[bool]=None) -> None:
        """Writes the coordinate cards in a sorted order"""
        model = self.model
        if model.coord.coord_id.max() > 0:
            bdf_file.write('$COORDS\n')
            model.coord.write_file(bdf_file, size=size, is_double=is_double)
        return

    def _write_dynamic(self, bdf_file: TextIOLike,
                       size: int=8, is_double: bool=False,
                       is_long_ids: Optional[bool]=None) -> None:
        """Writes the dynamic cards sorted by ID"""
        model = self.model
        model.tic.write_file(bdf_file, size=size, is_double=is_double, write_card_header=is_double)

        #bdf_file.write(model.tic.write(size=size))
        is_dynamic = (
            model.nlparms or model.frequencies or
            model.methods or model.cMethods or
            model.tsteps or model.tstepnls or
            #model.rotors
            model.nlpcis
        )
        if is_dynamic:
            bdf_file.write('$DYNAMIC\n')
            for (unused_id, method) in sorted(model.methods.items()):
                bdf_file.write(method.write_card(size, is_double))
            for (unused_id, cmethod) in sorted(model.cMethods.items()):
                bdf_file.write(cmethod.write_card(size, is_double))
            for (unused_id, nlparm) in sorted(model.nlparms.items()):
                bdf_file.write(nlparm.write_card(size, is_double))
            for (unused_id, nlpci) in sorted(model.nlpcis.items()):
                bdf_file.write(nlpci.write_card(size, is_double))
            for (unused_id, tstep) in sorted(model.tsteps.items()):
                bdf_file.write(tstep.write_card(size, is_double))
            for (unused_id, tstepnl) in sorted(model.tstepnls.items()):
                bdf_file.write(tstepnl.write_card(size, is_double))
            for (unused_id, freqs) in sorted(model.frequencies.items()):
                for freq in freqs:
                    bdf_file.write(freq.write_card(size, is_double))
            #for (unused_id, rotor) in sorted(model.rotors.items()):
                #bdf_file.write(rotor.write_card(size, is_double))

        #model.rotord.write_file(bdf_file, size=size, is_double=is_double)
        #model.rotor.write_file(bdf_file, size=size, is_double=is_double)
        model.rotorg.write_file(bdf_file, size=size, is_double=is_double)

    def _write_dmigs(self, bdf_file: TextIOLike,
                     size: int=8, is_double: bool=False,
                     is_long_ids: Optional[bool]=None) -> None:
        """
        Writes the DMIG cards

        Parameters
        ----------
        size : int
            large field (16) or small field (8)

        """
        model = self.model
        for (unused_name, dmig) in sorted(model.dmig.items()):
            bdf_file.write(dmig.write_card(size, is_double))
        for (unused_name, dmi) in sorted(model.dmi.items()):
            bdf_file.write(dmi.write_card(size, is_double))
        for (unused_name, dmij) in sorted(model.dmij.items()):
            bdf_file.write(dmij.write_card(size, is_double))
        for (unused_name, dmiji) in sorted(model.dmiji.items()):
            bdf_file.write(dmiji.write_card(size, is_double))
        for (unused_name, dmik) in sorted(model.dmik.items()):
            bdf_file.write(dmik.write_card(size, is_double))
        for (unused_name, dmiax) in sorted(model.dmiax.items()):
            bdf_file.write(dmiax.write_card(size, is_double))

    def _write_optimization(self, bdf_file: TextIOLike,
                            size: int=8, is_double: bool=False,
                            is_long_ids: Optional[bool]=None) -> None:
        """Writes the optimization cards sorted by ID"""
        model = self.model
        model.modtrak.write_file(bdf_file, size=size, is_double=is_double)
        model.ddval.write_file(bdf_file, size=size, is_double=is_double)
        model.desvar.write_file(bdf_file, size=size, is_double=is_double)
        model.dlink.write_file(bdf_file, size=size, is_double=is_double)
        model.dvgrid.write_file(bdf_file, size=size, is_double=is_double)

        model.dresp1.write_file(bdf_file, size=size, is_double=is_double)
        model.dresp2.write_file(bdf_file, size=size, is_double=is_double)  # poorly supported
        #model.dresp3.write_file(bdf_file, size=size, is_double=is_double)

        model.bndgrid.write_file(bdf_file, size=size, is_double=is_double)

        for deqatn_id, deqatn in model.dequations.items():
            bdf_file.write(deqatn.write_card(size=size, is_double=is_double))

        model.dconstr.write_file(bdf_file, size=size, is_double=is_double)
        model.dconadd.write_file(bdf_file, size=size, is_double=is_double)

        model.dvcrel1.write_file(bdf_file, size=size, is_double=is_double)
        model.dvcrel2.write_file(bdf_file, size=size, is_double=is_double)

        model.dvprel1.write_file(bdf_file, size=size, is_double=is_double)
        model.dvprel2.write_file(bdf_file, size=size, is_double=is_double)

        model.dvmrel1.write_file(bdf_file, size=size, is_double=is_double)
        model.dvmrel2.write_file(bdf_file, size=size, is_double=is_double)

        model.dscreen.write_file(bdf_file, size=size, is_double=is_double)
        if model.doptprm is not None:
            bdf_file.write(model.doptprm.write_card(size, is_double))
