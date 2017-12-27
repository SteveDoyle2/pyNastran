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
        self.beam = BEAM(self._h5n, self)
        self.bush = BUSH(self._h5n, self)
        self.quad4 = QUAD4(self._h5n, self)
        # TODO: what is better, quad4_cn or quad4cn?
        self.quad4_cn = QUAD4_CN(self._h5n, self)
        self.rod = ROD(self._h5n, self)
        self.shear = SHEAR(self._h5n, self)
        self.tria3 = TRIA3(self._h5n, self)

    def path(self):
        return self._elemental.path() + ['ELEMENT_FORCE']

########################################################################################################################


class BAR(ResultTable):
    """
    <dataset name="BAR">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="BM1A" type="double" description="Bending moment end A plane 1"/>
        <field name="BM2A" type="double" description="Bending moment end A plane 2"/>
        <field name="BM1B" type="double" description="Bending moment end B plane 1"/>
        <field name="BM2B" type="double" description="Bending moment end B plane 2"/>
        <field name="TS1" type="double" description="Shear plane 1"/>
        <field name="TS2" type="double" description="Shear plane 2"/>
        <field name="AF" type="double" description="Axial Force"/>
        <field name="TRQ" type="double" description="Torque"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 34 BAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR', result_type)

########################################################################################################################


class BEAM(ResultTable):
    """
    <dataset name="BEAM">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="FORCE" type="BEAM_FORCE" size="11" description="Element force structure for BEAM"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 2 BEAM REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM', result_type)

########################################################################################################################


class BUSH(ResultTable):
    """
    <dataset name="BUSH" sameAs="STRESS/BUSH"/>
    """

    result_type = 'ELEMENT FORCES 102 BUSH REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BUSH', result_type)

########################################################################################################################

class QUAD4(ResultTable):
    """
    <dataset name="QUAD4">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="FORCE" type="QUAD4_FORCE" description="data structure defined in typedef section"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 33 QUAD4 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4', result_type)


########################################################################################################################

class QUAD4_CN(ResultTable):
    """
    <dataset name="QUAD4_CN">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="TERM" type="character" size="4" description="Character string &quot;CEN/&quot;"/>
        <field name="FORCE" type="QUAD4_CN_FORCE" size="5" description="data structure defined in typedef section"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 144 QUAD4C REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN', result_type)

########################################################################################################################


class ROD(ResultTable):
    """
    <dataset name="ROD">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="AF" type="double" description="Axial Force"/>
        <field name="TRQ" type="double" description="Torque"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 1 ROD REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD', result_type)

########################################################################################################################


class SHEAR(ResultTable):
    """
    <dataset name="SHEAR">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="F41" type="double" description="Force  4 to 1"/>
        <field name="F21" type="double" description="Force  2 to 1"/>
        <field name="F12" type="double" description="Force  1 to 2"/>
        <field name="F32" type="double" description="Force  3 to 2"/>
        <field name="F23" type="double" description="Force  2 to 3"/>
        <field name="F43" type="double" description="Force  4 to 3"/>
        <field name="F34" type="double" description="Force  3 to 4"/>
        <field name="F14" type="double" description="Force  1 to 4"/>
        <field name="KF1" type="double" description="Kick Force on 1"/>
        <field name="S12" type="double" description="Shear  1    2"/>
        <field name="KF2" type="double" description="Kick Force on 2"/>
        <field name="S23" type="double" description="Shear  2    3"/>
        <field name="KF3" type="double" description="Kick Force on 3"/>
        <field name="S34" type="double" description="Shear  3    4"/>
        <field name="KF4" type="double" description="Kick Force on 4"/>
        <field name="S41" type="double" description="Shear  4    1"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 4 SHEAR REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/SHEAR', result_type)

########################################################################################################################


class TRIA3(ResultTable):
    """
    <dataset name="TRIA3">
        <field name="EID" type="integer" description="Element identification number"/>
        <field name="FORCE" type="QUAD4_FORCE" description="data structure defined in typedef section"/>
        <field name="DOMAIN_ID" type="integer" description="Domain identifier"/>
    </dataset>
    """

    result_type = 'ELEMENT FORCES 74 TRIA3 REAL'
    table_def = TableDef.create('/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3', result_type)
