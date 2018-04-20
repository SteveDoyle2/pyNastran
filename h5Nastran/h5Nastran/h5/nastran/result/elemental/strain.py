from __future__ import print_function, absolute_import

from h5Nastran.h5nastrannode import H5NastranNode
from ._shell_results import ShellElementStrainResultTable, ShellElementStrainResultTableComplex
from ..result_table import ResultTable, TableDef, DataGetter


class Strain(H5NastranNode):
    def __init__(self, h5n, elemental):
        self._h5n = h5n
        self._elemental = elemental

        self.axif2 = AXIF2(self._h5n, self)
        self.axif2_cplx = AXIF2_CPLX(self._h5n, self)
        self.axif3 = AXIF3(self._h5n, self)
        self.axif3_cplx = AXIF3_CPLX(self._h5n, self)
        self.axif4 = AXIF4(self._h5n, self)
        self.axif4_cplx = AXIF4_CPLX(self._h5n, self)
        # self.axisym = AXISYM(self._h5n, self)  # doesn't exist in MSC spec
        self.bar = BAR(self._h5n, self)
        self.bars = BARS(self._h5n, self)
        self.bars_cplx = BARS_CPLX(self._h5n, self)
        self.bar_cplx = BAR_CPLX(self._h5n, self)
        # self.bar_nl = BAR_NL(self._h5n, self)
        self.bar_rr = BAR_RR(self._h5n, self)
        self.beam = BEAM(self._h5n, self)
        self.beam3 = BEAM3(self._h5n, self)
        self.beam3_cplx = BEAM3_CPLX(self._h5n, self)
        self.beam_cplx = BEAM_CPLX(self._h5n, self)
        # self.beam_nl = BEAM_NL(self._h5n, self)
        self.beam_rr = BEAM_RR(self._h5n, self)
        self.bush = BUSH(self._h5n, self)
        self.bush1d = BUSH1D(self._h5n, self)
        self.bush1d_cplx = BUSH1D_CPLX(self._h5n, self)
        self.bush1d_rr = BUSH1D_RR(self._h5n, self)
        self.bush_cplx = BUSH_CPLX(self._h5n, self)
        # self.bush_nl = BUSH_NL(self._h5n, self)
        self.cone = CONE(self._h5n, self)
        self.conrod = CONROD(self._h5n, self)
        self.conrod_cplx = CONROD_CPLX(self._h5n, self)
        # self.conrod_nl = CONROD_NL(self._h5n, self)
        self.conrod_rr = CONROD_RR(self._h5n, self)
        self.elas1 = ELAS1(self._h5n, self)
        self.elas1_cplx = ELAS1_CPLX(self._h5n, self)
        # self.elas1_nl = ELAS1_NL(self._h5n, self)
        self.elas2 = ELAS2(self._h5n, self)
        self.elas2_cplx = ELAS2_CPLX(self._h5n, self)
        self.elas3 = ELAS3(self._h5n, self)
        self.elas3_cplx = ELAS3_CPLX(self._h5n, self)
        # self.elas3_nl = ELAS3_NL(self._h5n, self)
        # self.elem_comp = None  # skipping for now
        # self.extreme_fibre = None  # skipping for now
        # self.extreme_fibre_cplx = None  # skipping for now
        self.fast = FAST(self._h5n, self)
        self.fast_cplx = FAST_CPLX(self._h5n, self)
        self.gap = GAP(self._h5n, self)
        # self.gap_nl = GAP_NL(self._h5n, self)
        self.hexa = HEXA(self._h5n, self)
        # self.hexa20_27fdnl = HEXA20_27FDNL(self._h5n, self)
        # self.hexa20_8fdnl = HEXA20_8FDNL(self._h5n, self)
        # self.hexa20_fd = HEXA20_FD(self._h5n, self)
        self.hexa_cplx = HEXA_CPLX(self._h5n, self)
        # self.hexa_fd = HEXA_FD(self._h5n, self)
        # self.hexa_fdnl = HEXA_FDNL(self._h5n, self)
        # self.hexa_nl = HEXA_NL(self._h5n, self)
        # self.ifhexa = IFHEXA(self._h5n, self)
        # self.ifpenta = IFPENTA(self._h5n, self)
        self.penta = PENTA(self._h5n, self)
        # self.penta15_21fdnl = PENTA15_21FDNL(self._h5n, self)
        # self.penta15_6fdnl = PENTA15_6FDNL(self._h5n, self)
        # self.penta15_fd = PENTA15_FD(self._h5n, self)
        # self.penta_fdnl = PENTA_FDNL(self._h5n, self)
        # self.penta_nl = PENTA_NL(self._h5n, self)
        self.quad4 = QUAD4(self._h5n, self)
        self.quad4_comp = QUAD4_COMP(self._h5n, self)
        self.quad4_comp_cplx = QUAD4_COMP_CPLX(self._h5n, self)
        self.quad4_cplx = QUAD4_CPLX(self._h5n, self)
        self.quad4_fd = QUAD4_FD(self._h5n, self)
        # self.quad4_fdnl = QUAD4_FDNL(self._h5n, self)
        self.quad4_fd_cplx = QUAD4_FD_CPLX(self._h5n, self)
        # self.quad4_nl = QUAD4_NL(self._h5n, self)
        # self.quad8 = QUAD8(self._h5n, self)
        # self.quad8_4fdnl = QUAD8_4FDNL(self._h5n, self)
        # self.quad8_9fdnl = QUAD8_9FDNL(self._h5n, self)
        # self.quad8_comp = QUAD8_COMP(self._h5n, self)
        # self.quad8_comp_cplx = QUAD8_COMP_CPLX(self._h5n, self)
        # self.quad8_cplx = QUAD8_CPLX(self._h5n, self)
        self.quad8_fd = QUAD8_FD(self._h5n, self)
        self.quad8_fd_cplx = QUAD8_FD_CPLX(self._h5n, self)
        # self.quadr = QUADR(self._h5n, self)
        # self.quadr_comp = QUADR_COMP(self._h5n, self)
        # self.quadr_comp_cplx = QUADR_COMP_CPLX(self._h5n, self)
        # self.quadr_cplx = QUADR_CPLX(self._h5n, self)
        # self.quadr_fd = None  # skipping for now
        # self.quadr_fd_cplx = None  # skipping for now
        # self.quadr_nl = QUADR_NL(self._h5n, self)
        # self.quadx4_fd = QUADX4_FD(self._h5n, self)
        # self.quadx4_fdnl = QUADX4_FDNL(self._h5n, self)
        # self.quadx4_fd_cplx = QUADX4_FD_CPLX(self._h5n, self)
        # self.quadx8_4fdnl = QUADX8_4FDNL(self._h5n, self)
        # self.quadx8_9fdnl = QUADX8_9FDNL(self._h5n, self)
        # self.quadx8_fd = QUADX8_FD(self._h5n, self)
        # self.quadx8_fd_cplx = QUADX8_FD_CPLX(self._h5n, self)
        self.quad_cn = QUAD_CN(self._h5n, self)
        self.quad_cn_cplx = QUAD_CN_CPLX(self._h5n, self)
        # self.rac2d = RAC2D(self._h5n, self)
        # self.rac3d = RAC3D(self._h5n, self)
        self.rod = ROD(self._h5n, self)
        self.rod_cplx = ROD_CPLX(self._h5n, self)
        # self.rod_nl = ROD_NL(self._h5n, self)
        self.rod_rr = ROD_RR(self._h5n, self)
        # self.seam = SEAM(self._h5n, self)
        # self.seamp = SEAMP(self._h5n, self)
        # self.seamp_cplx = SEAMP_CPLX(self._h5n, self)
        self.shear = SHEAR(self._h5n, self)
        self.shear_cplx = SHEAR_CPLX(self._h5n, self)
        self.shear_rr = SHEAR_RR(self._h5n, self)
        # self.slot3 = SLOT3(self._h5n, self)
        # self.slot3_cplx = SLOT3_CPLX(self._h5n, self)
        # self.slot4 = SLOT4(self._h5n, self)
        # self.slot4_cplx = SLOT4_CPLX(self._h5n, self)
        self.tetra = TETRA(self._h5n, self)
        # self.tetra10_4fdnl = TETRA10_4FDNL(self._h5n, self)
        # self.tetra10_5fdnl = TETRA10_5FDNL(self._h5n, self)
        # self.tetra10_fd = TETRA10_FD(self._h5n, self)
        # self.tetra4_fdnl = TETRA4_FDNL(self._h5n, self)
        self.tetra_cplx = TETRA_CPLX(self._h5n, self)
        # self.tetra_fd = TETRA_FD(self._h5n, self)
        # self.tetra_fdnl = TETRA_FDNL(self._h5n, self)
        # self.tetra_nl = TETRA_NL(self._h5n, self)
        self.tria3 = TRIA3(self._h5n, self)
        # self.tria3_1fdnl = TRIA3_1FDNL(self._h5n, self)
        # self.tria3_3fdnl = TRIA3_3FDNL(self._h5n, self)
        # self.tria3_comp = TRIA3_COMP(self._h5n, self)
        # self.tria3_comp_cplx = TRIA3_COMP_CPLX(self._h5n, self)
        self.tria3_cplx = TRIA3_CPLX(self._h5n, self)
        self.tria3_fd = TRIA3_FD(self._h5n, self)
        self.tria3_fd_cplx = TRIA3_FD_CPLX(self._h5n, self)
        # self.tria3_nl = TRIA3_NL(self._h5n, self)
        self.tria6 = TRIA6(self._h5n, self)
        # self.tria6_comp = TRIA6_COMP(self._h5n, self)
        # self.tria6_comp_cplx = TRIA6_COMP_CPLX(self._h5n, self)
        self.tria6_cplx = TRIA6_CPLX(self._h5n, self)
        self.tria6_fd = TRIA6_FD(self._h5n, self)
        # self.tria6_fdnl = TRIA6_FDNL(self._h5n, self)
        self.tria6_fd_cplx = TRIA6_FD_CPLX(self._h5n, self)
        # self.triar = TRIAR(self._h5n, self)
        # self.triar_1fd = None  # skipping for now
        # self.triar_1fd_cplx = None  # skipping for now
        # self.triar_4fd = None  # skipping for now
        # self.triar_4fd_cplx = None  # skipping for now
        # self.triar_comp = TRIAR_COMP(self._h5n, self)
        # self.triar_comp_cplx = TRIAR_COMP_CPLX(self._h5n, self)
        # self.triar_cplx = TRIAR_CPLX(self._h5n, self)
        # self.triar_nl = TRIAR_NL(self._h5n, self)
        # self.triax3_1fdnl = TRIAX3_1FDNL(self._h5n, self)
        # self.triax3_3fdnl = TRIAX3_3FDNL(self._h5n, self)
        # self.triax3_fd = TRIAX3_FD(self._h5n, self)
        # self.triax3_fd_cplx = TRIAX3_FD_CPLX(self._h5n, self)
        self.triax6 = TRIAX6(self._h5n, self)
        self.triax6_cplx = TRIAX6_CPLX(self._h5n, self)
        # self.triax6_fd = TRIAX6_FD(self._h5n, self)
        # self.triax6_fdnl = TRIAX6_FDNL(self._h5n, self)
        # self.triax6_fd_cplx = TRIAX6_FD_CPLX(self._h5n, self)
        self.tube = TUBE(self._h5n, self)
        self.tube_cplx = TUBE_CPLX(self._h5n, self)
        # self.tube_nl = TUBE_NL(self._h5n, self)
        self.tube_rr = TUBE_RR(self._h5n, self)
        # self.visc_cplx = VISC_CPLX(self._h5n, self)
        # self.visc_rr = VISC_RR(self._h5n, self)
        self.weld = WELD(self._h5n, self)
        self.weld_cplx = WELD_CPLX(self._h5n, self)
        self.weldc = WELDC(self._h5n, self)
        self.weldc_cplx = WELDC_CPLX(self._h5n, self)
        self.weldp = WELDP(self._h5n, self)
        self.weldp_cplx = WELDP_CPLX(self._h5n, self)

    def path(self):
        return self._elemental.path() + ['STRAIN']

    def search(self, element_ids, domain_ids=()):
        # :type (List[int], List[int]) -> StrainResult
        result = StrainResult()
        _result = result.__dict__
        table_ids = self.__dict__.keys()
        _tables = self.__dict__
        for table_id in table_ids:
            if table_id.startswith('_'):  # not a table
                continue
            _result[table_id] = _tables[table_id].search(element_ids, domain_ids)

        return result

########################################################################################################################


class StrainResult(object):
    def __init__(self):
        self.axif2 = None  # type: DataFrame
        self.axif2_cplx = None  # type: DataFrame
        self.axif3 = None  # type: DataFrame
        self.axif3_cplx = None  # type: DataFrame
        self.axif4 = None  # type: DataFrame
        self.axif4_cplx = None  # type: DataFrame
        # self.axisym = None  # type: DataFrame
        self.bar = None  # type: DataFrame
        self.bars = None  # type: DataFrame
        self.bars_cplx = None  # type: DataFrame
        self.bar_cplx = None  # type: DataFrame
        # self.bar_nl = None  # type: DataFrame
        self.bar_rr = None  # type: DataFrame
        self.beam = None  # type: DataFrame
        self.beam3 = None  # type: DataFrame
        self.beam3_cplx = None  # type: DataFrame
        self.beam_cplx = None  # type: DataFrame
        # self.beam_nl = None  # type: DataFrame
        self.beam_rr = None  # type: DataFrame
        self.bush = None  # type: DataFrame
        self.bush1d = None  # type: DataFrame
        self.bush1d_cplx = None  # type: DataFrame
        self.bush1d_rr = None  # type: DataFrame
        self.bush_cplx = None  # type: DataFrame
        # self.bush_nl = None  # type: DataFrame
        self.cone = None  # type: DataFrame
        self.conrod = None  # type: DataFrame
        self.conrod_cplx = None  # type: DataFrame
        # self.conrod_nl = None  # type: DataFrame
        self.conrod_rr = None  # type: DataFrame
        self.elas1 = None  # type: DataFrame
        self.elas1_cplx = None  # type: DataFrame
        # self.elas1_nl = None  # type: DataFrame
        self.elas2 = None  # type: DataFrame
        self.elas2_cplx = None  # type: DataFrame
        self.elas3 = None  # type: DataFrame
        self.elas3_cplx = None  # type: DataFrame
        # self.elas3_nl = None  # type: DataFrame
        # self.elem_comp = None  # skipping for now
        # self.extreme_fibre = None  # skipping for now
        # self.extreme_fibre_cplx = None  # skipping for now
        self.fast = None  # type: DataFrame
        self.fast_cplx = None  # type: DataFrame
        self.gap = None  # type: DataFrame
        # self.gap_nl = None  # type: DataFrame
        self.hexa = None  # type: DataFrame
        # self.hexa20_27fdnl = None  # type: DataFrame
        # self.hexa20_8fdnl = None  # type: DataFrame
        # self.hexa20_fd = None  # type: DataFrame
        self.hexa_cplx = None  # type: DataFrame
        # self.hexa_fd = None  # type: DataFrame
        # self.hexa_fdnl = None  # type: DataFrame
        # self.hexa_nl = None  # type: DataFrame
        # self.ifhexa = None  # type: DataFrame
        # self.ifpenta = None  # type: DataFrame
        self.penta = None  # type: DataFrame
        # self.penta15_21fdnl = None  # type: DataFrame
        # self.penta15_6fdnl = None  # type: DataFrame
        # self.penta15_fd = None  # type: DataFrame
        # self.penta_fdnl = None  # type: DataFrame
        # self.penta_nl = None  # type: DataFrame
        self.quad4 = None  # type: DataFrame
        self.quad4_comp = None  # type: DataFrame
        self.quad4_comp_cplx = None  # type: DataFrame
        self.quad4_cplx = None  # type: DataFrame
        self.quad4_fd = None  # type: DataFrame
        # self.quad4_fdnl = None  # type: DataFrame
        self.quad4_fd_cplx = None  # type: DataFrame
        # self.quad4_nl = None  # type: DataFrame
        # self.quad8 = None  # type: DataFrame
        # self.quad8_4fdnl = None  # type: DataFrame
        # self.quad8_9fdnl = None  # type: DataFrame
        # self.quad8_comp = None  # type: DataFrame
        # self.quad8_comp_cplx = None  # type: DataFrame
        # self.quad8_cplx = None  # type: DataFrame
        self.quad8_fd = None  # type: DataFrame
        self.quad8_fd_cplx = None  # type: DataFrame
        # self.quadr = None  # type: DataFrame
        # self.quadr_comp = None  # type: DataFrame
        # self.quadr_comp_cplx = None  # type: DataFrame
        # self.quadr_cplx = None  # type: DataFrame
        # self.quadr_fd = None  # skipping for now
        # self.quadr_fd_cplx = None  # skipping for now
        # self.quadr_nl = None  # type: DataFrame
        # self.quadx4_fd = None  # type: DataFrame
        # self.quadx4_fdnl = None  # type: DataFrame
        # self.quadx4_fd_cplx = None  # type: DataFrame
        # self.quadx8_4fdnl = None  # type: DataFrame
        # self.quadx8_9fdnl = None  # type: DataFrame
        # self.quadx8_fd = None  # type: DataFrame
        # self.quadx8_fd_cplx = None  # type: DataFrame
        self.quad_cn = None  # type: DataFrame
        self.quad_cn_cplx = None  # type: DataFrame
        # self.rac2d = None  # type: DataFrame
        # self.rac3d = None  # type: DataFrame
        self.rod = None  # type: DataFrame
        self.rod_cplx = None  # type: DataFrame
        # self.rod_nl = None  # type: DataFrame
        self.rod_rr = None  # type: DataFrame
        # self.seam = None  # type: DataFrame
        # self.seamp = None  # type: DataFrame
        # self.seamp_cplx = None  # type: DataFrame
        self.shear = None  # type: DataFrame
        self.shear_cplx = None  # type: DataFrame
        self.shear_rr = None  # type: DataFrame
        # self.slot3 = None  # type: DataFrame
        # self.slot3_cplx = None  # type: DataFrame
        # self.slot4 = None  # type: DataFrame
        # self.slot4_cplx = None  # type: DataFrame
        self.tetra = None  # type: DataFrame
        # self.tetra10_4fdnl = None  # type: DataFrame
        # self.tetra10_5fdnl = None  # type: DataFrame
        # self.tetra10_fd = None  # type: DataFrame
        # self.tetra4_fdnl = None  # type: DataFrame
        self.tetra_cplx = None  # type: DataFrame
        # self.tetra_fd = None  # type: DataFrame
        # self.tetra_fdnl = None  # type: DataFrame
        # self.tetra_nl = None  # type: DataFrame
        self.tria3 = None  # type: DataFrame
        # self.tria3_1fdnl = None  # type: DataFrame
        # self.tria3_3fdnl = None  # type: DataFrame
        # self.tria3_comp = None  # type: DataFrame
        # self.tria3_comp_cplx = None  # type: DataFrame
        self.tria3_cplx = None  # type: DataFrame
        self.tria3_fd = None  # type: DataFrame
        self.tria3_fd_cplx = None  # type: DataFrame
        # self.tria3_nl = None  # type: DataFrame
        self.tria6 = None  # type: DataFrame
        # self.tria6_comp = None  # type: DataFrame
        # self.tria6_comp_cplx = None  # type: DataFrame
        self.tria6_cplx = None  # type: DataFrame
        self.tria6_fd = None  # type: DataFrame
        # self.tria6_fdnl = None  # type: DataFrame
        self.tria6_fd_cplx = None  # type: DataFrame
        # self.triar = None  # type: DataFrame
        # self.triar_1fd = None  # skipping for now
        # self.triar_1fd_cplx = None  # skipping for now
        # self.triar_4fd = None  # skipping for now
        # self.triar_4fd_cplx = None  # skipping for now
        # self.triar_comp = None  # type: DataFrame
        # self.triar_comp_cplx = None  # type: DataFrame
        # self.triar_cplx = None  # type: DataFrame
        # self.triar_nl = None  # type: DataFrame
        # self.triax3_1fdnl = None  # type: DataFrame
        # self.triax3_3fdnl = None  # type: DataFrame
        # self.triax3_fd = None  # type: DataFrame
        # self.triax3_fd_cplx = None  # type: DataFrame
        self.triax6 = None  # type: DataFrame
        self.triax6_cplx = None  # type: DataFrame
        # self.triax6_fd = None  # type: DataFrame
        # self.triax6_fdnl = None  # type: DataFrame
        # self.triax6_fd_cplx = None  # type: DataFrame
        self.tube = None  # type: DataFrame
        self.tube_cplx = None  # type: DataFrame
        # self.tube_nl = None  # type: DataFrame
        self.tube_rr = None  # type: DataFrame
        # self.visc_cplx = None  # type: DataFrame
        # self.visc_rr = None  # type: DataFrame
        self.weld = None  # type: DataFrame
        self.weld_cplx = None  # type: DataFrame
        self.weldc = None  # type: DataFrame
        self.weldc_cplx = None  # type: DataFrame
        self.weldp = None  # type: DataFrame
        self.weldp_cplx = None  # type: DataFrame


########################################################################################################################


class AXIF2(ResultTable):
    result_type = 'ELEMENT STRAINS 47 AXIF2 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF2', result_type)


########################################################################################################################


class AXIF2_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 47 AXIF2 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF2_CPLX', result_type)


########################################################################################################################


class AXIF3(ResultTable):
    result_type = 'ELEMENT STRAINS 48 AXIF3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF3', result_type)


########################################################################################################################


class AXIF3_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 48 AXIF3 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF3_CPLX', result_type)


########################################################################################################################


class AXIF4(ResultTable):
    result_type = 'ELEMENT STRAINS 49 AXIF4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF4', result_type)


########################################################################################################################


class AXIF4_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 49 AXIF4 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXIF4_CPLX', result_type)


########################################################################################################################

# doesn't exist in MSC spec
# class AXISYM(ResultTable):
#     result_type = 'ELEMENT STRAINS 241 AXISYM REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/AXISYM', result_type)


########################################################################################################################


class BAR(ResultTable):
    result_type = 'ELEMENT STRAINS 34 BAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BAR', result_type)


########################################################################################################################


class BARS(ResultTable):
    result_type = 'ELEMENT STRAINS 100 BARS REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BARS', result_type)


########################################################################################################################


class BARS_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 100 BARS COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BARS_CPLX', result_type)


########################################################################################################################


class BAR_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 34 BAR COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BAR_CPLX', result_type)


########################################################################################################################


# class BAR_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 240 BAR_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BAR_NL', result_type)


########################################################################################################################


class BAR_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 34 BAR_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BAR_RR', result_type)


########################################################################################################################


class BEAM(ResultTable):
    result_type = 'ELEMENT STRAINS 2 BEAM REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM', result_type)


########################################################################################################################


class BEAM3(ResultTable):
    result_type = 'ELEMENT STRAINS 184 BEAM3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM3', result_type)


########################################################################################################################


class BEAM3_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 184 BEAM3 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM3_CPLX', result_type)


########################################################################################################################


class BEAM_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 2 BEAM COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM_CPLX', result_type)


########################################################################################################################


# class BEAM_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 94 BEAM_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM_NL', result_type)


########################################################################################################################


class BEAM_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 2 BEAM_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BEAM_RR', result_type)


########################################################################################################################


class BUSH(ResultTable):
    result_type = 'ELEMENT STRAINS 102 BUSH REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH', result_type)


########################################################################################################################


class BUSH1D(ResultTable):
    result_type = 'ELEMENT STRAINS 40 BUSH1D REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH1D', result_type)


########################################################################################################################


class BUSH1D_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 40 BUSH1D COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH1D_CPLX', result_type)


########################################################################################################################


class BUSH1D_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 40 BUSH1D_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH1D_RR', result_type)


########################################################################################################################


class BUSH_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 102 BUSH COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH_CPLX', result_type)


########################################################################################################################


# class BUSH_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 226 BUSH_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/BUSH_NL', result_type)


########################################################################################################################


class CONE(ResultTable):
    result_type = 'ELEMENT STRAINS 35 CONE REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/CONE', result_type)


########################################################################################################################


class CONROD(ResultTable):
    result_type = 'ELEMENT STRAINS 10 CONROD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/CONROD', result_type)


########################################################################################################################


class CONROD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 10 CONROD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/CONROD_CPLX', result_type)


########################################################################################################################


# class CONROD_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 92 CONROD_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/CONROD_NL', result_type)


########################################################################################################################


class CONROD_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 10 CONROD_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/CONROD_RR', result_type)


########################################################################################################################


class ELAS1(ResultTable):
    result_type = 'ELEMENT STRAINS 11 ELAS1 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS1', result_type)


########################################################################################################################


class ELAS1_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 11 ELAS1 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS1_CPLX', result_type)


########################################################################################################################


# class ELAS1_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 224 ELAS1_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS1_NL', result_type)


########################################################################################################################


class ELAS2(ResultTable):
    result_type = 'ELEMENT STRAINS 12 ELAS2 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS2', result_type)


########################################################################################################################


class ELAS2_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 12 ELAS2 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS2_CPLX', result_type)


########################################################################################################################


class ELAS3(ResultTable):
    result_type = 'ELEMENT STRAINS 13 ELAS3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS3', result_type)


########################################################################################################################


class ELAS3_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 13 ELAS3 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS3_CPLX', result_type)


########################################################################################################################


# class ELAS3_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 225 ELAS3_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELAS3_NL', result_type)


########################################################################################################################


class FAST(ResultTable):
    result_type = 'ELEMENT STRAINS 126 FAST REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/FAST', result_type)


########################################################################################################################


class FAST_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 126 FAST COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/FAST_CPLX', result_type)


########################################################################################################################


class GAP(ResultTable):
    result_type = 'ELEMENT STRAINS 38 GAP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/GAP', result_type)


########################################################################################################################


# class GAP_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 86 GAP_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/GAP_NL', result_type)


########################################################################################################################


class HEXA(ResultTable):
    result_type = 'ELEMENT STRAINS 67 HEXA REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA', result_type)


########################################################################################################################


# class HEXA20_27FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 207 HEXA20_27FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA20_27FDNL', result_type)


########################################################################################################################


# class HEXA20_8FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 202 HEXA20_8FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA20_8FDNL', result_type)


########################################################################################################################


# class HEXA20_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 163 HEXA20_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA20_FD', result_type)


########################################################################################################################


class HEXA_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 67 HEXA COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA_CPLX', result_type)


########################################################################################################################


# class HEXA_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 140 HEXA_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA_FD', result_type)


########################################################################################################################


# class HEXA_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 218 HEXA_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA_FDNL', result_type)


########################################################################################################################


# class HEXA_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 93 HEXA_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/HEXA_NL', result_type)


########################################################################################################################


# class IFHEXA(ResultTable):
#     result_type = 'ELEMENT STRAINS 65 IFHEXA REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/IFHEXA', result_type)


########################################################################################################################


# class IFPENTA(ResultTable):
#     result_type = 'ELEMENT STRAINS 66 IFPENTA REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/IFPENTA', result_type)


########################################################################################################################


class PENTA(ResultTable):
    result_type = 'ELEMENT STRAINS 68 PENTA REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA', result_type)


########################################################################################################################


# class PENTA15_21FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 209 PENTA15_21FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA15_21FDNL', result_type)


########################################################################################################################


# class PENTA15_6FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 204 PENTA15_6FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA15_6FDNL', result_type)


########################################################################################################################


# class PENTA15_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 165 PENTA15_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA15_FD', result_type)


########################################################################################################################


# class PENTA_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 220 PENTA_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA_FDNL', result_type)


########################################################################################################################


# class PENTA_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 91 PENTA_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/PENTA_NL', result_type)


########################################################################################################################


class QUAD4(ResultTable, ShellElementStrainResultTable):
    result_type = 'ELEMENT STRAINS 33 QUAD4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4', result_type)
    table_def.add_index_option('MATERIAL', None)
    table_def.add_index_option('FIBER', None)
    table_def.add_index_option('VONM', DataGetter(indices=[0, 2, 3, 4, 5, 10, 11, 12, 13]))


########################################################################################################################


class QUAD4_COMP(ResultTable):
    result_type = 'ELEMENT STRAINS 95 QUAD4LC REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_COMP', result_type)


########################################################################################################################


class QUAD4_COMP_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 95 QUAD4LC COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_COMP_CPLX', result_type)


########################################################################################################################


class QUAD4_CPLX(ResultTable, ShellElementStrainResultTableComplex):
    result_type = 'ELEMENT STRAINS 33 QUAD4 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_CPLX', result_type)


########################################################################################################################


class QUAD4_FD(ResultTable):
    result_type = 'ELEMENT STRAINS 139 QUAD4_FD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_FD', result_type)


########################################################################################################################


# class QUAD4_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 201 QUAD4_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_FDNL', result_type)


########################################################################################################################


class QUAD4_FD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 139 QUAD4_FD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_FD_CPLX', result_type)


########################################################################################################################


# class QUAD4_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 90 QUAD4_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD4_NL', result_type)


########################################################################################################################


# class QUAD8(ResultTable):
#     result_type = 'ELEMENT STRAINS 64 QUAD8 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8', result_type)


########################################################################################################################


# class QUAD8_4FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 219 QUAD8_4FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_4FDNL', result_type)


########################################################################################################################


# class QUAD8_9FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 208 QUAD8_9FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_9FDNL', result_type)


########################################################################################################################


# class QUAD8_COMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 96 QUAD8LC REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_COMP', result_type)


########################################################################################################################


# class QUAD8_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 96 QUAD8LC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_COMP_CPLX', result_type)


########################################################################################################################


# class QUAD8_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 64 QUAD8 COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_CPLX', result_type)


########################################################################################################################


class QUAD8_FD(ResultTable):
    result_type = 'ELEMENT STRAINS 164 QUAD8_FD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_FD', result_type)


########################################################################################################################


class QUAD8_FD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 164 QUAD8_FD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD8_FD_CPLX', result_type)


########################################################################################################################


# class QUADR(ResultTable):
#     result_type = 'ELEMENT STRAINS 82 QUADR REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADR', result_type)


########################################################################################################################


# class QUADR_COMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 232 QUADRLC REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADR_COMP', result_type)


########################################################################################################################


# class QUADR_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 232 QUADRLC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADR_COMP_CPLX', result_type)


########################################################################################################################


# class QUADR_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 82 QUADR COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADR_CPLX', result_type)


########################################################################################################################


# class QUADR_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 172 QUADR_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADR_NL', result_type)


########################################################################################################################


# class QUADX4_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 170 QUADX4_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX4_FD', result_type)


########################################################################################################################


# class QUADX4_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 214 QUADX4_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX4_FDNL', result_type)


########################################################################################################################


# class QUADX4_FD_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 170 QUADX4_FD COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX4_FD_CPLX', result_type)


########################################################################################################################


# class QUADX8_4FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 223 QUADX8_4FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX8_4FDNL', result_type)


########################################################################################################################


# class QUADX8_9FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 215 QUADX8_9FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX8_9FDNL', result_type)


########################################################################################################################


# class QUADX8_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 171 QUADX8_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX8_FD', result_type)


########################################################################################################################


# class QUADX8_FD_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 171 QUADX8_FD COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUADX8_FD_CPLX', result_type)


########################################################################################################################


class QUAD_CN(ResultTable):
    result_type = 'ELEMENT STRAINS 144 QUADC REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD_CN', result_type)


########################################################################################################################


class QUAD_CN_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 144 QUADC COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/QUAD_CN_CPLX', result_type)


########################################################################################################################


# class RAC2D(ResultTable):
#     result_type = 'ELEMENT STRAINS 60 RAC2D REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/RAC2D', result_type)


########################################################################################################################


# class RAC3D(ResultTable):
#     result_type = 'ELEMENT STRAINS 61 RAC3D REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/RAC3D', result_type)


########################################################################################################################


class ROD(ResultTable):
    result_type = 'ELEMENT STRAINS 1 ROD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ROD', result_type)


########################################################################################################################


class ROD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 1 ROD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ROD_CPLX', result_type)


########################################################################################################################


# class ROD_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 89 ROD_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ROD_NL', result_type)


########################################################################################################################


class ROD_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 1 ROD_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/ROD_RR', result_type)


########################################################################################################################


# class SEAM(ResultTable):
#     result_type = 'ELEMENT STRAINS 119 SEAM REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SEAM', result_type)


########################################################################################################################


# class SEAMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 159 SEAMP REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SEAMP', result_type)


########################################################################################################################


# class SEAMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 159 SEAMP COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SEAMP_CPLX', result_type)


########################################################################################################################


class SHEAR(ResultTable):
    result_type = 'ELEMENT STRAINS 4 SHEAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SHEAR', result_type)


########################################################################################################################


class SHEAR_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 4 SHEAR COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SHEAR_CPLX', result_type)


########################################################################################################################


class SHEAR_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 4 SHEAR_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SHEAR_RR', result_type)


########################################################################################################################


# class SLOT3(ResultTable):
#     result_type = 'ELEMENT STRAINS 50 SLOT3 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SLOT3', result_type)


########################################################################################################################


# class SLOT3_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 50 SLOT3 COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SLOT3_CPLX', result_type)


########################################################################################################################


# class SLOT4(ResultTable):
#     result_type = 'ELEMENT STRAINS 51 SLOT4 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SLOT4', result_type)


########################################################################################################################


# class SLOT4_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 51 SLOT4 COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/SLOT4_CPLX', result_type)


########################################################################################################################


class TETRA(ResultTable):
    result_type = 'ELEMENT STRAINS 39 TETRA REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA', result_type)


########################################################################################################################


# class TETRA10_4FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 221 TETRA10_4FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA10_4FDNL', result_type)


########################################################################################################################


# class TETRA10_5FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 210 TETRA10_5FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA10_5FDNL', result_type)


########################################################################################################################


# class TETRA10_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 166 TETRA10_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA10_FD', result_type)


########################################################################################################################


# class TETRA4_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 216 TETRA4_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA4_FDNL', result_type)


########################################################################################################################


class TETRA_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 39 TETRA COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA_CPLX', result_type)


########################################################################################################################


# class TETRA_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 161 TETRA_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA_FD', result_type)


########################################################################################################################


# class TETRA_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 205 TETRA_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA_FDNL', result_type)


########################################################################################################################


# class TETRA_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 85 TETRA_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TETRA_NL', result_type)


########################################################################################################################


class TRIA3(ResultTable, ShellElementStrainResultTable):
    result_type = 'ELEMENT STRAINS 74 TRIA3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3', result_type)
    table_def.add_index_option('MATERIAL', None)
    table_def.add_index_option('FIBER', None)
    table_def.add_index_option('VONM', DataGetter(indices=[0, 2, 3, 4, 5, 10, 11, 12, 13]))

########################################################################################################################


# class TRIA3_1FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 206 TRIA3_1FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_1FDNL', result_type)


########################################################################################################################


# class TRIA3_3FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 217 TRIA3_3FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_3FDNL', result_type)


########################################################################################################################


# class TRIA3_COMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 97 TRIA3LC REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_COMP', result_type)


########################################################################################################################


# class TRIA3_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 97 TRIA3LC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_COMP_CPLX', result_type)


########################################################################################################################


class TRIA3_CPLX(ResultTable, ShellElementStrainResultTableComplex):
    result_type = 'ELEMENT STRAINS 74 TRIA3 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_CPLX', result_type)


########################################################################################################################


class TRIA3_FD(ResultTable):
    result_type = 'ELEMENT STRAINS 162 TRIA3_FD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_FD', result_type)


########################################################################################################################


class TRIA3_FD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 162 TRIA3_FD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_FD_CPLX', result_type)


########################################################################################################################


# class TRIA3_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 88 TRIA3_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA3_NL', result_type)


########################################################################################################################


class TRIA6(ResultTable, ShellElementStrainResultTable):
    result_type = 'ELEMENT STRAINS 75 TRIA6 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6', result_type)


########################################################################################################################


# class TRIA6_COMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 98 TRIA6LC REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_COMP', result_type)


########################################################################################################################


# class TRIA6_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 98 TRIA6LC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_COMP_CPLX', result_type)


########################################################################################################################


class TRIA6_CPLX(ResultTable, ShellElementStrainResultTableComplex):
    result_type = 'ELEMENT STRAINS 75 TRIA6 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_CPLX', result_type)


########################################################################################################################


class TRIA6_FD(ResultTable):
    result_type = 'ELEMENT STRAINS 167 TRIA6_FD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_FD', result_type)


########################################################################################################################


# class TRIA6_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 211 TRIA6_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_FDNL', result_type)


########################################################################################################################


class TRIA6_FD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 167 TRIA6_FD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIA6_FD_CPLX', result_type)


########################################################################################################################


# class TRIAR(ResultTable):
#     result_type = 'ELEMENT STRAINS 70 TRIAR REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAR', result_type)


########################################################################################################################


# class TRIAR_COMP(ResultTable):
#     result_type = 'ELEMENT STRAINS 233 TRIARLC REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAR_COMP', result_type)


########################################################################################################################


# class TRIAR_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 233 TRIARLC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAR_COMP_CPLX', result_type)


########################################################################################################################


# class TRIAR_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 70 TRIAR COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAR_CPLX', result_type)


########################################################################################################################


# class TRIAR_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 173 TRIAR_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAR_NL', result_type)


########################################################################################################################


# class TRIAX3_1FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 212 TRIAX3_1FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX3_1FDNL', result_type)


########################################################################################################################


# class TRIAX3_3FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 222 TRIAX3_3FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX3_3FDNL', result_type)


########################################################################################################################


# class TRIAX3_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 168 TRIAX3_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX3_FD', result_type)


########################################################################################################################


# class TRIAX3_FD_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 168 TRIAX3_FD COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX3_FD_CPLX', result_type)


########################################################################################################################


class TRIAX6(ResultTable):
    result_type = 'ELEMENT STRAINS 53 TRIAX6 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX6', result_type)


########################################################################################################################


class TRIAX6_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 53 TRIAX6 COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX6_CPLX', result_type)


########################################################################################################################


# class TRIAX6_FD(ResultTable):
#     result_type = 'ELEMENT STRAINS 169 TRIAX6_FD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX6_FD', result_type)


########################################################################################################################


# class TRIAX6_FDNL(ResultTable):
#     result_type = 'ELEMENT STRAINS 213 TRIAX6_FDNL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX6_FDNL', result_type)


########################################################################################################################


# class TRIAX6_FD_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 169 TRIAX6_FD COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TRIAX6_FD_CPLX', result_type)


########################################################################################################################


class TUBE(ResultTable):
    result_type = 'ELEMENT STRAINS 3 TUBE REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TUBE', result_type)


########################################################################################################################


class TUBE_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 3 TUBE COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TUBE_CPLX', result_type)


########################################################################################################################


# class TUBE_NL(ResultTable):
#     result_type = 'ELEMENT STRAINS 87 TUBE_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TUBE_NL', result_type)


########################################################################################################################


class TUBE_RR(ResultTable):
    result_type = 'ELEMENT STRAINS 3 TUBE_RR RANDOM'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/TUBE_RR', result_type)


########################################################################################################################


# class VISC_CPLX(ResultTable):
#     result_type = 'ELEMENT STRAINS 24 VISC COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/VISC_CPLX', result_type)


########################################################################################################################


# class VISC_RR(ResultTable):
#     result_type = 'ELEMENT STRAINS 24 VISC_RR RANDOM'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/VISC_RR', result_type)


########################################################################################################################


class WELD(ResultTable):
    result_type = 'ELEMENT STRAINS 200 WELD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELD', result_type)


########################################################################################################################


class WELD_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 200 WELD COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELD_CPLX', result_type)


########################################################################################################################


class WELDC(ResultTable):
    result_type = 'ELEMENT STRAINS 117 WELDC REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELDC', result_type)


########################################################################################################################


class WELDC_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 117 WELDC COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELDC_CPLX', result_type)


########################################################################################################################


class WELDP(ResultTable):
    result_type = 'ELEMENT STRAINS 118 WELDP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELDP', result_type)


########################################################################################################################


class WELDP_CPLX(ResultTable):
    result_type = 'ELEMENT STRAINS 118 WELDP COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/STRAIN/WELDP_CPLX', result_type)


