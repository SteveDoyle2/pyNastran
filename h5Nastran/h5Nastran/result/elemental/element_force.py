from __future__ import print_function, absolute_import
from six import iteritems, itervalues
from six.moves import range

import tables
import numpy as np

from ..result_table import ResultTable, TableDef


class ElementForce(object):
    def __init__(self, h5n, elemental):
        self._h5n = h5n
        self._elemental = elemental

        self.bar = BAR(self._h5n, self)
        self.bar_cplx = BAR_CPLX(self._h5n, self)
        self.bars = BARS(self._h5n, self)
        self.bars_cplx = BARS_CPLX(self._h5n, self)
        self.beam = BEAM(self._h5n, self)
        self.beam_cplx = BEAM_CPLX(self._h5n, self)
        self.beam3 = BEAM3(self._h5n, self)
        self.beam3_cplx = BEAM3_CPLX(self._h5n, self)
        self.bend = BEND(self._h5n, self)
        self.bend_cplx = BEND_CPLX(self._h5n, self)
        self.bush = BUSH(self._h5n, self)
        self.bush_cplx = BUSH_CPLX(self._h5n, self)
        # self.coneax = CONEAX(self._h5n, self)
        # self.coneax_cplx = CONEAX_CPLX(self._h5n, self)
        self.damp1 = DAMP1(self._h5n, self)
        self.damp1_cplx = DAMP1_CPLX(self._h5n, self)
        self.damp2 = DAMP2(self._h5n, self)
        self.damp2_cplx = DAMP2_CPLX(self._h5n, self)
        self.damp3 = DAMP3(self._h5n, self)
        self.damp3_cplx = DAMP3_CPLX(self._h5n, self)
        self.damp4 = DAMP4(self._h5n, self)
        self.damp4_cplx = DAMP4_CPLX(self._h5n, self)
        # self.dum3 = DUM3(self._h5n, self)
        # self.dum3_cplx = DUM3_CPLX(self._h5n, self)
        # self.dum4 = DUM4(self._h5n, self)
        # self.dum4_cplx = DUM4_CPLX(self._h5n, self)
        # self.dum5 = DUM5(self._h5n, self)
        # self.dum5_cplx = DUM5_CPLX(self._h5n, self)
        # self.dum6 = DUM6(self._h5n, self)
        # self.dum6_cplx = DUM6_CPLX(self._h5n, self)
        # self.dum7 = DUM7(self._h5n, self)
        # self.dum7_cplx = DUM7_CPLX(self._h5n, self)
        # self.dum8 = DUM8(self._h5n, self)
        # self.dum8_cplx = DUM8_CPLX(self._h5n, self)
        # self.dum9 = DUM9(self._h5n, self)
        # self.dum9_cplx = DUM9_CPLX(self._h5n, self)
        self.elas1 = ELAS1(self._h5n, self)
        self.elas1_cplx = ELAS1_CPLX(self._h5n, self)
        self.elas2 = ELAS2(self._h5n, self)
        self.elas2_cplx = ELAS2_CPLX(self._h5n, self)
        self.elas3 = ELAS3(self._h5n, self)
        self.elas3_cplx = ELAS3_CPLX(self._h5n, self)
        self.elas4 = ELAS4(self._h5n, self)
        self.elas4_cplx = ELAS4_CPLX(self._h5n, self)
        self.gap = GAP(self._h5n, self)
        # self.gap_cplx = GAP_CPLX(self._h5n, self)
        self.conrod = CONROD(self._h5n, self)
        self.conrod_cplx = CONROD_CPLX(self._h5n, self)
        self.quad4 = QUAD4(self._h5n, self)
        self.quad4_cplx = QUAD4_CPLX(self._h5n, self)
        self.quad4_comp = QUAD4_COMP(self._h5n, self)
        # self.quad4_comp_cplx = QUAD4_COMP_CPLX(self._h5n, self)
        self.quad4_cn = QUAD4_CN(self._h5n, self)
        self.quad4_cn_cplx = QUAD4_CN_CPLX(self._h5n, self)
        self.quad8 = QUAD8(self._h5n, self)
        self.quad8_cplx = QUAD8_CPLX(self._h5n, self)
        self.quad8_comp = QUAD8_COMP(self._h5n, self)
        # self.quad8_comp_cplx = QUAD8_COMP_CPLX(self._h5n, self)
        self.quadr = QUADR(self._h5n, self)
        self.quadr_cplx = QUADR_CPLX(self._h5n, self)
        # self.quadr_nl = QUADR_NL(self._h5n, self)
        # self.quadr_nl_cplx = QUADR_NL_CPLX(self._h5n, self)
        self.rod = ROD(self._h5n, self)
        self.rod_cplx = ROD_CPLX(self._h5n, self)
        self.shear = SHEAR(self._h5n, self)
        self.shear_cplx = SHEAR_CPLX(self._h5n, self)
        self.tria3 = TRIA3(self._h5n, self)
        self.tria3_cplx = TRIA3_CPLX(self._h5n, self)
        self.tria3_comp = TRIA3_COMP(self._h5n, self)
        # self.tria3_comp_cplx = TRIA3_COMP_CPLX(self._h5n, self)
        self.tria6 = TRIA6(self._h5n, self)
        self.tria6_cplx = TRIA6_CPLX(self._h5n, self)
        self.tria6_comp = TRIA6_COMP(self._h5n, self)
        # self.tria6_comp_cplx = TRIA6_COMP_CPLX(self._h5n, self)
        self.triar = TRIAR(self._h5n, self)
        self.triar_cplx = TRIAR_CPLX(self._h5n, self)
        # self.triar_nl = TRIAR_NL(self._h5n, self)
        # self.triar_nl_cplx = TRIAR_NL_CPLX(self._h5n, self)
        self.tube = TUBE(self._h5n, self)
        self.tube_cplx = TUBE_CPLX(self._h5n, self)
        self.visc = VISC(self._h5n, self)
        self.visc_cplx = VISC_CPLX(self._h5n, self)
        self.weldp = WELDP(self._h5n, self)
        self.weldp_cplx = WELDP_CPLX(self._h5n, self)
        self.weldc = WELDC(self._h5n, self)
        self.weldc_cplx = WELDC_CPLX(self._h5n, self)
        self.weld = WELD(self._h5n, self)
        self.weld_cplx = WELD_CPLX(self._h5n, self)
        # self.vuquad = VUQUAD(self._h5n, self)
        # self.vuquad_cplx = VUQUAD_CPLX(self._h5n, self)
        # self.vutria = VUTRIA(self._h5n, self)
        # self.vutria_cplx = VUTRIA_CPLX(self._h5n, self)
        # self.vubeam = VUBEAM(self._h5n, self)
        # self.vubeam_cplx = VUBEAM_CPLX(self._h5n, self)

    def path(self):
        return self._elemental.path() + ['ELEMENT_FORCE']


########################################################################################################################


class BAR(ResultTable):
    result_type = 'ELEMENT FORCES 34 BAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR', result_type)


########################################################################################################################


class BAR_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 34 BAR_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX', result_type)


########################################################################################################################


class BARS(ResultTable):
    result_type = 'ELEMENT FORCES 100 BARS REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BARS', result_type)


########################################################################################################################


class BARS_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 100 BARS_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BARS_CPLX', result_type)


########################################################################################################################


class BEAM(ResultTable):
    result_type = 'ELEMENT FORCES 2 BEAM REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM', result_type)


########################################################################################################################


class BEAM_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 2 BEAM_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM_CPLX', result_type)


########################################################################################################################


class BEAM3(ResultTable):
    result_type = 'ELEMENT FORCES 184 BEAM3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM3', result_type)


########################################################################################################################


class BEAM3_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 184 BEAM3_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM3_CPLX', result_type)


########################################################################################################################


class BEND(ResultTable):
    result_type = 'ELEMENT FORCES 69 BEND REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEND', result_type)


########################################################################################################################


class BEND_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 69 BEND_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEND_CPLX', result_type)


########################################################################################################################


class BUSH(ResultTable):
    result_type = 'ELEMENT FORCES 102 BUSH REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BUSH', result_type)


########################################################################################################################


class BUSH_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 102 BUSH_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BUSH_CPLX', result_type)


########################################################################################################################


# class CONEAX(ResultTable):
#     result_type = 'ELEMENT FORCES 35 CONEAX REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/CONEAX', result_type)
#
#
# ########################################################################################################################
#
#
# class CONEAX_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 35 CONEAX_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/CONEAX_CPLX', result_type)


########################################################################################################################


class DAMP1(ResultTable):
    result_type = 'ELEMENT FORCES 20 DAMP1 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP1', result_type)


########################################################################################################################


class DAMP1_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 20 DAMP1_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP1_CPLX', result_type)


########################################################################################################################


class DAMP2(ResultTable):
    result_type = 'ELEMENT FORCES 21 DAMP2 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP2', result_type)


########################################################################################################################


class DAMP2_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 21 DAMP2_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP2_CPLX', result_type)


########################################################################################################################


class DAMP3(ResultTable):
    result_type = 'ELEMENT FORCES 22 DAMP3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP3', result_type)


########################################################################################################################


class DAMP3_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 22 DAMP3_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP3_CPLX', result_type)


########################################################################################################################


class DAMP4(ResultTable):
    result_type = 'ELEMENT FORCES 23 DAMP4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP4', result_type)


########################################################################################################################


class DAMP4_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 23 DAMP4_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP4_CPLX', result_type)


########################################################################################################################


# class DUM3(ResultTable):
#     result_type = 'ELEMENT FORCES 55 DUM3 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM3', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM3_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 55 DUM3_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM3_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM4(ResultTable):
#     result_type = 'ELEMENT FORCES 56 DUM4 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM4', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM4_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 56 DUM4_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM4_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM5(ResultTable):
#     result_type = 'ELEMENT FORCES 57 DUM5 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM5', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM5_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 57 DUM5_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM5_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM6(ResultTable):
#     result_type = 'ELEMENT FORCES 58 DUM6 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM6', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM6_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 58 DUM6_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM6_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM7(ResultTable):
#     result_type = 'ELEMENT FORCES 59 DUM7 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM7', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM7_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 59 DUM7_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM7_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM8(ResultTable):
#     result_type = 'ELEMENT FORCES 60 DUM8 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM8', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM8_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 60 DUM8_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM8_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM9(ResultTable):
#     result_type = 'ELEMENT FORCES 61 DUM9 REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM9', result_type)
#
#
# ########################################################################################################################
#
#
# class DUM9_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 61 DUM9_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/DUM9_CPLX', result_type)
#

########################################################################################################################


class ELAS1(ResultTable):
    result_type = 'ELEMENT FORCES 11 ELAS1 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1', result_type)


########################################################################################################################


class ELAS1_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 11 ELAS1_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX', result_type)


########################################################################################################################


class ELAS2(ResultTable):
    result_type = 'ELEMENT FORCES 12 ELAS2 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS2', result_type)


########################################################################################################################


class ELAS2_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 12 ELAS2_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS2_CPLX', result_type)


########################################################################################################################


class ELAS3(ResultTable):
    result_type = 'ELEMENT FORCES 13 ELAS3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS3', result_type)


########################################################################################################################


class ELAS3_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 13 ELAS3_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS3_CPLX', result_type)


########################################################################################################################


class ELAS4(ResultTable):
    result_type = 'ELEMENT FORCES 14 ELAS4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS4', result_type)


########################################################################################################################


class ELAS4_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 14 ELAS4_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS4_CPLX', result_type)


########################################################################################################################


class GAP(ResultTable):
    result_type = 'ELEMENT FORCES 38 GAP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/GAP', result_type)


########################################################################################################################


# class GAP_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 38 GAP_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/GAP_CPLX', result_type)


########################################################################################################################


class CONROD(ResultTable):
    result_type = 'ELEMENT FORCES 10 CONROD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/CONROD', result_type)


########################################################################################################################


class CONROD_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 10 CONROD_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/CONROD_CPLX', result_type)


########################################################################################################################


class QUAD4(ResultTable):
    result_type = 'ELEMENT FORCES 33 QUAD4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4', result_type)


########################################################################################################################


class QUAD4_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 33 QUAD4_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CPLX', result_type)


########################################################################################################################


class QUAD4_COMP(ResultTable):
    result_type = 'ELEMENT FORCES 95 QUAD4_COMP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_COMP', result_type)


########################################################################################################################


# class QUAD4_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 95 QUAD4_COMP_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_COMP_CPLX', result_type)


########################################################################################################################


class QUAD4_CN(ResultTable):
    result_type = 'ELEMENT FORCES 144 QUAD4C REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN', result_type)


########################################################################################################################


class QUAD4_CN_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 144 QUAD4C COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN_CPLX', result_type)


########################################################################################################################


class QUAD8(ResultTable):
    result_type = 'ELEMENT FORCES 64 QUAD8 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8', result_type)


########################################################################################################################


class QUAD8_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 64 QUAD8_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8_CPLX', result_type)


########################################################################################################################


class QUAD8_COMP(ResultTable):
    result_type = 'ELEMENT FORCES 96 QUAD8_COMP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8_COMP', result_type)


########################################################################################################################


# class QUAD8_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 96 QUAD8_COMP_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8_COMP_CPLX', result_type)


########################################################################################################################


class QUADR(ResultTable):
    result_type = 'ELEMENT FORCES 82 QUADR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR', result_type)


########################################################################################################################


class QUADR_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 82 QUADR_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR_CPLX', result_type)


########################################################################################################################


# class QUADR_NL(ResultTable):
#     result_type = 'ELEMENT FORCES 235 QUADR_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR_NL', result_type)


########################################################################################################################


# class QUADR_NL_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 235 QUADR_NL_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR_NL_CPLX', result_type)


########################################################################################################################


class ROD(ResultTable):
    result_type = 'ELEMENT FORCES 1 ROD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD', result_type)


########################################################################################################################


class ROD_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 1 ROD_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD_CPLX', result_type)


########################################################################################################################


class SHEAR(ResultTable):
    result_type = 'ELEMENT FORCES 4 SHEAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/SHEAR', result_type)


########################################################################################################################


class SHEAR_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 4 SHEAR_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/SHEAR_CPLX', result_type)


########################################################################################################################


class TRIA3(ResultTable):
    result_type = 'ELEMENT FORCES 74 TRIA3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3', result_type)


########################################################################################################################


class TRIA3_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 74 TRIA3_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3_CPLX', result_type)


########################################################################################################################


class TRIA3_COMP(ResultTable):
    result_type = 'ELEMENT FORCES 97 TRIA3_COMP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3_COMP', result_type)


########################################################################################################################


# class TRIA3_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 97 TRIA3_COMP_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3_COMP_CPLX', result_type)


########################################################################################################################


class TRIA6(ResultTable):
    result_type = 'ELEMENT FORCES 75 TRIA6 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6', result_type)


########################################################################################################################


class TRIA6_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 75 TRIA6_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6_CPLX', result_type)


########################################################################################################################


class TRIA6_COMP(ResultTable):
    result_type = 'ELEMENT FORCES 98 TRIA6_COMP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6_COMP', result_type)


########################################################################################################################


# class TRIA6_COMP_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 98 TRIA6_COMP_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6_COMP_CPLX', result_type)


########################################################################################################################


class TRIAR(ResultTable):
    result_type = 'ELEMENT FORCES 70 TRIAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR', result_type)


########################################################################################################################


class TRIAR_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 70 TRIAR_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR_CPLX', result_type)


########################################################################################################################


# class TRIAR_NL(ResultTable):
#     result_type = 'ELEMENT FORCES 236 TRIAR_NL REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR_NL', result_type)


########################################################################################################################


# class TRIAR_NL_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 236 TRIAR_NL_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR_NL_CPLX', result_type)


########################################################################################################################


class TUBE(ResultTable):
    result_type = 'ELEMENT FORCES 3 TUBE REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TUBE', result_type)


########################################################################################################################


class TUBE_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 3 TUBE_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TUBE_CPLX', result_type)


########################################################################################################################


class VISC(ResultTable):
    result_type = 'ELEMENT FORCES 24 VISC REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VISC', result_type)


########################################################################################################################


class VISC_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 24 VISC_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VISC_CPLX', result_type)


########################################################################################################################


class WELDP(ResultTable):
    result_type = 'ELEMENT FORCES 118 WELDP REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELDP', result_type)


########################################################################################################################


class WELDP_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 118 WELDP_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELDP_CPLX', result_type)


########################################################################################################################


class WELDC(ResultTable):
    result_type = 'ELEMENT FORCES 117 WELDC REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELDC', result_type)


########################################################################################################################


class WELDC_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 117 WELDC_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELDC_CPLX', result_type)


########################################################################################################################


class WELD(ResultTable):
    result_type = 'ELEMENT FORCES 200 WELD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELD', result_type)


########################################################################################################################


class WELD_CPLX(ResultTable):
    result_type = 'ELEMENT FORCES 200 WELD_CPLX COMPLEX'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/WELD_CPLX', result_type)


########################################################################################################################


# class VUQUAD(ResultTable):
#     result_type = 'ELEMENT FORCES 189 VUQUAD REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUQUAD', result_type)
#
#
# ########################################################################################################################
#
#
# class VUQUAD_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 189 VUQUAD_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUQUAD_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class VUTRIA(ResultTable):
#     result_type = 'ELEMENT FORCES 190 VUTRIA REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUTRIA', result_type)
#
#
# ########################################################################################################################
#
#
# class VUTRIA_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 190 VUTRIA_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUTRIA_CPLX', result_type)
#
#
# ########################################################################################################################
#
#
# class VUBEAM(ResultTable):
#     result_type = 'ELEMENT FORCES 191 VUBEAM REAL'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUBEAM', result_type)
#
#
# ########################################################################################################################
#
#
# class VUBEAM_CPLX(ResultTable):
#     result_type = 'ELEMENT FORCES 191 VUBEAM_CPLX COMPLEX'
#     table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/VUBEAM_CPLX', result_type)


