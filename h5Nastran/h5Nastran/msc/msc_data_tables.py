from collections import OrderedDict

data_tables = OrderedDict()


def register_table(table):
    table_id = table.path + "/" + table.name
    data_tables[table_id] = table
    return table


@register_table
class TERMS(object):
    name = 'TERMS'
    path = '/NASTRAN/INPUT/CONSTRAINT/AELINK'
    dtype = [('LABLI', 'S8', ()), ('CI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/AELINK'
    dtype = [('ID', '<i8', ()), ('LABLD', 'S8', ()), ('TERMS_POS', '<i8', ()), ('TERMS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/AELINK/TERMS']


@register_table
class CSSCHD(object):
    name = 'CSSCHD'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('AESID', '<i8', ()), ('LALPHA', '<i8', ()), ('LMACH', '<i8', ()),
             ('LSCHD', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CYSUP(object):
    name = 'CYSUP'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('GID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DEFORM(object):
    name = 'DEFORM'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRDSET(object):
    name = 'GRDSET'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('ID', '<i8', ()), ('CP', '<i8', ()), ('CD', '<i8', ()), ('PS', '<i8', ()), ('SEID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GCA(object):
    name = 'GCA'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPC'
    dtype = [('G', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPC'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ()), ('GCA_POS', '<i8', ()),
             ('GCA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/MPC/GCA']


@register_table
class S(object):
    name = 'S'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCADD'
    dtype = [('S', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCADD'
    dtype = [('SID', '<i8', ()), ('S_POS', '<i8', ()), ('S_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/MPCADD/S']


@register_table
class TERMS(object):
    name = 'TERMS'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCAX'
    dtype = [('RID1', '<i8', ()), ('HID1', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCAX'
    dtype = [('SID', '<i8', ()), ('RID', '<i8', ()), ('HID', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ()),
             ('TERMS_POS', '<i8', ()), ('TERMS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/MPCAX/TERMS']


@register_table
class MPCD(object):
    name = 'MPCD'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('GM', '<i8', ()), ('CM', '<i8', ()), ('Y', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCY'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('AI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/MPCY'
    dtype = [('SID', '<i8', ()), ('GM', '<i8', ()), ('CM', '<i8', ()), ('AM', '<f8', ()), ('Y', '<f8', ()),
             ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/MPCY/GRIDS']


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/CONSTRAINT/RSPLINE'
    dtype = [('G2', '<i8', ()), ('C2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/RSPLINE'
    dtype = [('EID', '<i8', ()), ('DBYL', '<f8', ()), ('G1', '<i8', ()), ('G_POS', '<i8', ()), ('G_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/RSPLINE/G']


@register_table
class SESUP(object):
    name = 'SESUP'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SEID', '<i8', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPBLND1(object):
    name = 'SPBLND1'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('SID1', '<i8', ()), ('SID2', '<i8', ()), ('OPT', 'S8', ()), ('W1', '<f8', ()),
             ('GID', '<i8', ()), ('D1', '<f8', ()), ('D2', '<f8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()),
             ('X3', '<f8', ()), ('CID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPBLND2(object):
    name = 'SPBLND2'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('SID1', '<i8', ()), ('SID2', '<i8', ()), ('OPT', 'S8', ()), ('AELIST', '<i8', ()),
             ('D1', '<f8', ()), ('D2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPC(object):
    name = 'SPC'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_G'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_G'
    dtype = [('SID', '<i8', ()), ('C', '<i8', ()), ('G_POS', '<i8', ()), ('G_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SPC1/SPC1_G/G']


@register_table
class SPC1_THRU(object):
    name = 'SPC1_THRU'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPC1'
    dtype = [('SID', '<i8', ()), ('C', '<i8', ()), ('FIRST', '<i8', ()), ('SECOND', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class S(object):
    name = 'S'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPCADD'
    dtype = [('S', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPCADD'
    dtype = [('SID', '<i8', ()), ('S_POS', '<i8', ()), ('S_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SPCADD/S']


@register_table
class SPCAX(object):
    name = 'SPCAX'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('RID', '<i8', ()), ('HID', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPCD(object):
    name = 'SPCD'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPCOFF(object):
    name = 'SPCOFF'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('G', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GIS(object):
    name = 'GIS'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPCOFF1'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPCOFF1'
    dtype = [('C', '<i8', ()), ('G_POS', '<i8', ()), ('G_LEN', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SPCOFF1/GIS']


@register_table
class SPCR(object):
    name = 'SPCR'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINE1(object):
    name = 'SPLINE1'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('BOX1', '<i8', ()), ('BOX2', '<i8', ()), ('SETG', '<i8', ()),
             ('DZ', '<f8', ()), ('METHOD', 'S8', ()), ('USAGE', 'S8', ()), ('NELEM', '<i8', ()), ('MELEM', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINE2(object):
    name = 'SPLINE2'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('SETG', '<i8', ()),
             ('DZ', '<f8', ()), ('DTOR', '<f8', ()), ('CID', '<i8', ()), ('DTHX', '<f8', ()), ('DTHY', '<f8', ()),
             ('USAGE', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GCA(object):
    name = 'GCA'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPLINE3'
    dtype = [('G', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPLINE3'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('UKID', '<i8', ()), ('COMP', '<i8', ()), ('USAGE', 'S8', ()),
             ('GCA_POS', '<i8', ()), ('GCA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SPLINE3/GCA']


@register_table
class SPLINE4(object):
    name = 'SPLINE4'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('AELIST', '<i8', ()), ('SETG', '<i8', ()), ('DZ', '<f8', ()),
             ('METHOD', 'S8', ()), ('USAGE', 'S8', ()), ('NELEM', '<i8', ()), ('MELEM', '<i8', ()),
             ('FTYPE', '<i8', ()), ('RCORE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINE5(object):
    name = 'SPLINE5'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('AELIST', '<i8', ()), ('SETG', '<i8', ()), ('DZ', '<f8', ()),
             ('DTORXY', '<f8', ()), ('CID', '<i8', ()), ('DTHX', '<f8', ()), ('DTHY', '<f8', ()), ('DTHZ', '<f8', ()),
             ('USAGE', 'S8', ()), ('METHOD', 'S8', ()), ('DTORZY', '<f8', ()), ('FTYPE', '<i8', ()),
             ('RCORE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINE6(object):
    name = 'SPLINE6'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('AELIST', '<i8', ()), ('SETG', '<i8', ()), ('DZ', '<f8', ()),
             ('METHOD', 'S8', ()), ('USAGE', 'S8', ()), ('VSTYPE', 'S8', ()), ('VSLIST', '<i8', ()),
             ('I2VNUM', '<i8', ()), ('D2VNUM', '<i8', ()), ('METHVS', 'S8', ()), ('DZR', '<f8', ()),
             ('METHCON', 'S8', ()), ('NGRID', '<i8', ()), ('ELTOL', '<f8', ()), ('NCYCLE', '<i8', ()),
             ('AUGWEI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINE7(object):
    name = 'SPLINE7'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('CAERO', '<i8', ()), ('AELIST', '<i8', ()), ('SETG', '<i8', ()), ('DZ', '<f8', ()),
             ('DTOR', '<f8', ()), ('CID', '<i8', ()), ('USAGE', 'S8', ()), ('METHOD', 'S8', ()), ('DZR', '<f8', ()),
             ('IA2', '<f8', ()), ('EPSBM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPLINEX(object):
    name = 'SPLINEX'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('EID', '<i8', ()), ('GROUP', 'S8', ()), ('DGCOMP', 'S8', ()), ('IGCOMP', 'S8', ()), ('DECOMP', 'S8', ()),
             ('IECOMP', 'S8', ()), ('USAGE', 'S8', ()), ('AELIST', '<i8', ()), ('AEFACT', '<i8', ()),
             ('AELISTC', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPLINRB'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SPLINRB'
    dtype = [('SID', '<i8', ()), ('CAERO', '<i8', ()), ('AELIST', '<i8', ()), ('USAGE', 'S8', ()),
             ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SPLINRB/GRIDS']


@register_table
class SUPAX(object):
    name = 'SUPAX'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('RID', '<i8', ()), ('HID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SUPORT(object):
    name = 'SUPORT'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COMPONENT(object):
    name = 'COMPONENT'
    path = '/NASTRAN/INPUT/CONSTRAINT/SUPORT1'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/SUPORT1'
    dtype = [('SID', '<i8', ()), ('COMPONENT_POS', '<i8', ()), ('COMPONENT_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/SUPORT1/COMPONENT']


@register_table
class TEMP(object):
    name = 'TEMP'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPAX(object):
    name = 'TEMPAX'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('RID', '<i8', ()), ('PHI', '<f8', ()), ('TEMP', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPB3(object):
    name = 'TEMPB3'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('TA', '<f8', ()), ('TB', '<f8', ()), ('TC', '<f8', ()),
             ('TPYA', '<f8', ()), ('TPZA', '<f8', ()), ('TPYB', '<f8', ()), ('TPZB', '<f8', ()), ('TPYC', '<f8', ()),
             ('TPZC', '<f8', ()), ('TCA', '<f8', ()), ('TDA', '<f8', ()), ('TEA', '<f8', ()), ('TFA', '<f8', ()),
             ('TCB', '<f8', ()), ('TDB', '<f8', ()), ('TEB', '<f8', ()), ('TFB', '<f8', ()), ('TCC', '<f8', ()),
             ('TDC', '<f8', ()), ('TEC', '<f8', ()), ('TFC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPBC(object):
    name = 'TEMPBC'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S4', ()), ('TEMP1', '<f8', ()), ('GID1', '<i8', ()), ('TEMP2', '<f8', ()),
             ('GID2', '<i8', ()), ('TEMP3', '<f8', ()), ('GID3', '<i8', ()), ('GIDII', '<i8', ()), ('INC', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPD(object):
    name = 'TEMPD'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPN1(object):
    name = 'TEMPN1'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPP1(object):
    name = 'TEMPP1'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('TBAR', '<f8', ()), ('TPRIME', '<f8', ()), ('T1', '<f8', ()),
             ('T2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPP2(object):
    name = 'TEMPP2'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('T', '<f8', ()), ('MX', '<f8', ()), ('MY', '<f8', ()),
             ('MXY', '<f8', ()), ('TT', '<f8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPP3(object):
    name = 'TEMPP3'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('Z', '<f8', (11,)), ('T', '<f8', (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPRB(object):
    name = 'TEMPRB'
    path = '/NASTRAN/INPUT/CONSTRAINT'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('TA', '<f8', ()), ('TB', '<f8', ()), ('TP1A', '<f8', ()),
             ('TP1B', '<f8', ()), ('TP2A', '<f8', ()), ('TP2B', '<f8', ()), ('TS', '<f8', (8,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIMS(object):
    name = 'TRIMS'
    path = '/NASTRAN/INPUT/CONSTRAINT/TRIM'
    dtype = [('LABEL', 'S8', ()), ('UX', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/TRIM'
    dtype = [('ID', '<i8', ()), ('MACH', '<f8', ()), ('Q', '<f8', ()), ('AEQR', '<f8', ()), ('TRIMS_POS', '<i8', ()),
             ('TRIMS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/TRIM/TRIMS']


@register_table
class STATES(object):
    name = 'STATES'
    path = '/NASTRAN/INPUT/CONSTRAINT/TRIM2'
    dtype = [('LABEL', 'S8', ()), ('CLASS', 'S8', ()), ('UX', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/TRIM2'
    dtype = [('ID', '<i8', ()), ('MACH', '<f8', ()), ('Q', '<f8', ()), ('AEQR', '<f8', ()), ('STATES_POS', '<i8', ()),
             ('STATES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/TRIM2/STATES']


@register_table
class STATES(object):
    name = 'STATES'
    path = '/NASTRAN/INPUT/CONSTRAINT/UXVEC'
    dtype = [('LABEL', 'S8', ()), ('UX', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONSTRAINT/UXVEC'
    dtype = [('ID', '<i8', ()), ('STATES_POS', '<i8', ()), ('STATES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONSTRAINT/UXVEC/STATES']


@register_table
class BCBDPRP(object):
    name = 'BCBDPRP'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('ID', '<i8', ()), ('BNCI', '<i8', ()), ('BNCR', '<f8', ()), ('BNCEI', '<i8', ()), ('BNCER', '<f8', ()),
             ('BNLI', '<i8', ()), ('BNLR', '<f8', ()), ('BNLEI', '<i8', ()), ('BNLER', '<f8', ()),
             ('CFILMI', '<i8', ()), ('CFILMR', '<f8', ()), ('CMB', '<f8', ()), ('CMS', '<f8', ()), ('COPTB', '<i8', ()),
             ('EMISSI', '<i8', ()), ('EMISSR', '<f8', ()), ('FRICI', '<i8', ()), ('FRICR', '<f8', ()),
             ('HBLI', '<i8', ()), ('HBLR', '<f8', ()), ('HCTI', '<i8', ()), ('HCTR', '<f8', ()), ('HCVI', '<i8', ()),
             ('HCVR', '<f8', ()), ('HNCI', '<i8', ()), ('HNCR', '<f8', ()), ('HNCEI', '<i8', ()), ('HNCER', '<f8', ()),
             ('HNLI', '<i8', ()), ('HNLR', '<f8', ()), ('HNLEI', '<i8', ()), ('HNLER', '<f8', ()), ('IDSPL', '<i8', ()),
             ('ISTYP', '<i8', ()), ('ITYPE', '<i8', ()), ('MIDNOD', '<i8', ()), ('SANGLE', '<f8', ()),
             ('TBODYI', '<i8', ()), ('TBODYR', '<f8', ()), ('TSINKI', '<i8', ()), ('TSINKR', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/CONTACT/BCBMRAD'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/CONTACT/BCBMRAD'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU_BY(object):
    name = 'THRU_BY'
    path = '/NASTRAN/INPUT/CONTACT/BCBMRAD'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('N', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCBMRAD'
    dtype = [('RADIUS', '<f8', ()), ('TYPE', 'S4', ()), ('ALL', '<i8', ()), ('LIST_POS', '<i8', ()),
             ('LIST_LEN', '<i8', ()), ('THRU_POS', '<i8', ()), ('THRU_LEN', '<i8', ()), ('THRUBY_POS', '<i8', ()),
             ('THRUBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCBMRAD/IDLIST', '/NASTRAN/INPUT/CONTACT/BCBMRAD/THRU',
                 '/NASTRAN/INPUT/CONTACT/BCBMRAD/THRU_BY']


@register_table
class ADVANCE(object):
    name = 'ADVANCE'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('SANGEL', '<f8', ()), ('COPTB', '<i8', ()), ('USER', 'S4', ()),
             ('MIDNODE', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class APPROV(object):
    name = 'APPROV'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('AI', '<i8', ()), ('AR', '<f8', ()), ('N1I', '<i8', ()), ('N1R', '<f8', ()),
             ('N2I', '<i8', ()), ('N2R', '<f8', ()), ('N3I', '<i8', ()), ('N3R', '<f8', ()), ('V1I', '<i8', ()),
             ('V1R', '<f8', ()), ('V2I', '<i8', ()), ('V2R', '<f8', ()), ('V3I', '<i8', ()), ('V3R', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class ARC(object):
    name = 'ARC'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPTSA', '<i8', ()), ('METHOD', '<i8', ()), ('GRID_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class BEZIER(object):
    name = 'BEZIER'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NP1', '<i8', ()), ('NP2', '<i8', ()), ('NSUB1', '<i8', ()), ('NSUB2', '<i8', ()),
             ('GRIDS_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class COORD2D(object):
    name = 'COORD2D'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class COORDS(object):
    name = 'COORDS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class CPOINTS(object):
    name = 'CPOINTS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('CPOINT', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class CYLIND(object):
    name = 'CYLIND'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NSUB', '<i8', ()), ('GTOP', '<i8', ()), ('RTOP', '<f8', ()), ('GBOTTOM', '<i8', ()),
             ('RBOTTOM', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GENERAL(object):
    name = 'GENERAL'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('NLOAD', '<i8', ()), ('ANGVELI', '<i8', ()), ('ANGVELR', '<f8', ()), ('DCOS1I', '<i8', ()),
             ('DCOS1R', '<f8', ()), ('DCOS2I', '<i8', ()), ('DCOS2R', '<f8', ()), ('DCOS3I', '<i8', ()),
             ('DCOS3R', '<f8', ()), ('VELRB1I', '<i8', ()), ('VELRB1R', '<f8', ()), ('VELRB2I', '<i8', ()),
             ('VELRB2R', '<f8', ()), ('VELRB3I', '<i8', ()), ('VELRB3R', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('GRID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GROW(object):
    name = 'GROW'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('GF1', '<f8', ()), ('GF2', '<f8', ()), ('GF3', '<f8', ()), ('TABGF1', '<i8', ()),
             ('TABGF2', '<i8', ()), ('TABGF3', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class HEAT(object):
    name = 'HEAT'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('CFILMI', '<i8', ()), ('CFILMR', '<f8', ()), ('TSINKI', '<i8', ()),
             ('TSINKR', '<f8', ()), ('CHEATI', '<i8', ()), ('CHEATR', '<f8', ()), ('TBODYI', '<i8', ()),
             ('TBODYR', '<f8', ()), ('HCVI', '<i8', ()), ('HCVR', '<f8', ()), ('HNCI', '<i8', ()), ('HNCR', '<f8', ()),
             ('ITYPE', '<i8', ()), ('BNCI', '<i8', ()), ('BNCR', '<f8', ()), ('EMISSI', '<i8', ()),
             ('EMISSR', '<f8', ()), ('HBLI', '<i8', ()), ('HBLR', '<f8', ()), ('HNLI', '<i8', ()), ('HNLR', '<f8', ()),
             ('BNLI', '<i8', ()), ('BNLR', '<f8', ()), ('HNLEI', '<i8', ()), ('HNLER', '<f8', ()), ('BNLEI', '<i8', ()),
             ('BNLER', '<f8', ()), ('HNCEI', '<i8', ()), ('HNCER', '<f8', ()), ('BNCEI', '<i8', ()),
             ('BNCER', '<f8', ()), ('CMBI', '<i8', ()), ('CMBR', '<f8', ()), ('CMSI', '<i8', ()), ('CMSR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class HOMOCOORS(object):
    name = 'HOMOCOORS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('HOMOCOOR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class KNOTS(object):
    name = 'KNOTS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('KNOT', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class LINE(object):
    name = 'LINE'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPTSL', '<i8', ()), ('GRID_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class NURBS(object):
    name = 'NURBS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPTU', '<i8', ()), ('NPTV', '<i8', ()), ('NORU', '<i8', ()), ('NORV', '<i8', ()),
             ('NSUBU', '<i8', ()), ('NSUBV', '<i8', ()), ('NTRIM', '<i8', ()), ('CPOINT_FLAG', '<i8', ()),
             ('CPOINT_POS', '<i8', ()), ('COORD_POS', '<i8', ()), ('HOMOCOOR_POS', '<i8', ()), ('KNOT_POS', '<i8', ()),
             ('TRIMVEC_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class NURBS2(object):
    name = 'NURBS2'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('IDN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class NURBS2D(object):
    name = 'NURBS2D'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPTU2D', '<i8', ()), ('NORU2D', '<i8', ()), ('NSUB', '<i8', ()),
             ('COORD_FLAG', '<i8', ()), ('CPOINT_POS', '<i8', ()), ('COORD2D_POS', '<i8', ()),
             ('HOMOCOORD_POS', '<i8', ()), ('KNOT_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PATCH3D(object):
    name = 'PATCH3D'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPATCH', '<i8', ()), ('PATCHES_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PATCHES(object):
    name = 'PATCHES'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('IDP', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()), ('G4', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class POLY(object):
    name = 'POLY'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NP1', '<i8', ()), ('NP2', '<i8', ()), ('GRIDS_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class RIGID(object):
    name = 'RIGID'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('CGID', '<i8', ()), ('NENT', '<i8', ()), ('NAME', 'S40', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SPHERE(object):
    name = 'SPHERE'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NSUB', '<i8', ()), ('GCENTER', '<i8', ()), ('RADIUS', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SPLINE(object):
    name = 'SPLINE'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('TYPE', 'S8', ()), ('NPTSS', '<i8', ()), ('GRID_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRIMHOMOS(object):
    name = 'TRIMHOMOS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('HOMOCOOR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRIMISOS(object):
    name = 'TRIMISOS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('XISO', '<f8', ()), ('YISO', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRIMKNOTS(object):
    name = 'TRIMKNOTS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('KNOT', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRIMVECS(object):
    name = 'TRIMVECS'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('IDTRIM', '<i8', ()), ('NPTUTRIM', '<i8', ()), ('NORUTRIM', '<i8', ()), ('NSUBTRIM', '<i8', ()),
             ('TRIMISO_POS', '<i8', ()), ('TRIMHOMO_POS', '<i8', ()), ('TRIMKNOT_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCBODY'
    dtype = [('BID', '<i8', ()), ('DIM', '<i8', ()), ('BEHAV', 'S8', ()), ('BSID', '<i8', ()), ('ISTYP', '<i8', ()),
             ('FRICI', '<i8', ()), ('FRICR', '<f8', ()), ('IDSPL', '<i8', ()), ('CONTROL', '<i8', ()),
             ('GENERAL_POS', '<i8', ()), ('GENERAL_LEN', '<i8', ()), ('RIGID_POS', '<i8', ()), ('RIGID_LEN', '<i8', ()),
             ('HEAT_POS', '<i8', ()), ('HEAT_LEN', '<i8', ()), ('PATCH3D_POS', '<i8', ()), ('PATCH3D_LEN', '<i8', ()),
             ('BEZIER_POS', '<i8', ()), ('BEZIER_LEN', '<i8', ()), ('POLY_POS', '<i8', ()), ('POLY_LEN', '<i8', ()),
             ('CYLIND_POS', '<i8', ()), ('CYLIND_LEN', '<i8', ()), ('SPHERE_POS', '<i8', ()), ('SPHERE_LEN', '<i8', ()),
             ('NURBS2_POS', '<i8', ()), ('NURBS2_LEN', '<i8', ()), ('NURBS_POS', '<i8', ()), ('NURBS_LEN', '<i8', ()),
             ('NURBS2D_POS', '<i8', ()), ('NURBS2D_LEN', '<i8', ()), ('LINE_POS', '<i8', ()), ('LINE_LEN', '<i8', ()),
             ('ARC_POS', '<i8', ()), ('ARC_LEN', '<i8', ()), ('SPLINE_POS', '<i8', ()), ('SPLINE_LEN', '<i8', ()),
             ('ADVANCE_POS', '<i8', ()), ('ADVANCE_LEN', '<i8', ()), ('APPROV_POS', '<i8', ()),
             ('APPROV_LEN', '<i8', ()), ('GROW_POS', '<i8', ()), ('GROW_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCBODY/ADVANCE', '/NASTRAN/INPUT/CONTACT/BCBODY/APPROV',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/ARC', '/NASTRAN/INPUT/CONTACT/BCBODY/BEZIER',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/COORD2D', '/NASTRAN/INPUT/CONTACT/BCBODY/COORDS',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/CPOINTS', '/NASTRAN/INPUT/CONTACT/BCBODY/CYLIND',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/GENERAL', '/NASTRAN/INPUT/CONTACT/BCBODY/GRIDS',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/GROW', '/NASTRAN/INPUT/CONTACT/BCBODY/HEAT',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/HOMOCOORS', '/NASTRAN/INPUT/CONTACT/BCBODY/KNOTS',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/LINE', '/NASTRAN/INPUT/CONTACT/BCBODY/NURBS',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/NURBS2', '/NASTRAN/INPUT/CONTACT/BCBODY/NURBS2D',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/PATCH3D', '/NASTRAN/INPUT/CONTACT/BCBODY/PATCHES',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/POLY', '/NASTRAN/INPUT/CONTACT/BCBODY/RIGID',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/SPHERE', '/NASTRAN/INPUT/CONTACT/BCBODY/SPLINE',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/TRIMHOMOS', '/NASTRAN/INPUT/CONTACT/BCBODY/TRIMISOS',
                 '/NASTRAN/INPUT/CONTACT/BCBODY/TRIMKNOTS', '/NASTRAN/INPUT/CONTACT/BCBODY/TRIMVECS']


@register_table
class BCBODY1(object):
    name = 'BCBODY1'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('BID', '<i8', ()), ('BPID', '<i8', ()), ('DIM', 'S8', ()), ('BEHAV', 'S8', ()), ('BSID', '<i8', ()),
             ('BCGRID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/CONTACT/BCBZIER'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCBZIER'
    dtype = [('RBID', '<i8', ()), ('NP1', '<i8', ()), ('NP2', '<i8', ()), ('NSUB1', '<i8', ()), ('NSUB2', '<i8', ()),
             ('NGRID', '<i8', ()), ('FIRST', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCBZIER/G']


@register_table
class RIGIDBODY(object):
    name = 'RIGIDBODY'
    path = '/NASTRAN/INPUT/CONTACT/BCHANGE'
    dtype = [('IDRBOD', '<i8', ()), ('N1', '<i8', ()), ('N2', '<i8', ()), ('INC', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCHANGE'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('NBOD', '<i8', ()), ('RB_POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCHANGE/RIGIDBODY']


@register_table
class RIGIDBODY(object):
    name = 'RIGIDBODY'
    path = '/NASTRAN/INPUT/CONTACT/BCMOVE'
    dtype = [('IDRBOD', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCMOVE'
    dtype = [('ID', '<i8', ()), ('MTYPE', 'S8', ()), ('IREL', '<i8', ()), ('RBOD_POS', '<i8', ()),
             ('RBOD_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCMOVE/RIGIDBODY']


@register_table
class GRID(object):
    name = 'GRID'
    path = '/NASTRAN/INPUT/CONTACT/BCNURB2'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class HOMO(object):
    name = 'HOMO'
    path = '/NASTRAN/INPUT/CONTACT/BCNURB2'
    dtype = [('HOMOI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class KNOT(object):
    name = 'KNOT'
    path = '/NASTRAN/INPUT/CONTACT/BCNURB2'
    dtype = [('KNOTI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/CONTACT/BCNURB2'
    dtype = [('XI', '<f8', ()), ('YI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCNURB2'
    dtype = [('RBID', '<i8', ()), ('NPTU', '<i8', ()), ('NORU', '<i8', ()), ('NSUB', '<i8', ()), ('NGRID', '<i8', ()),
             ('NXY', '<i8', ()), ('NHOMO', '<i8', ()), ('NKNOT', '<i8', ()), ('GRID_POS', '<i8', ()),
             ('XY_POS', '<i8', ()), ('HOMO_POS', '<i8', ()), ('KNOT_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCNURB2/GRID', '/NASTRAN/INPUT/CONTACT/BCNURB2/HOMO',
                 '/NASTRAN/INPUT/CONTACT/BCNURB2/KNOT', '/NASTRAN/INPUT/CONTACT/BCNURB2/XY']


@register_table
class GRID(object):
    name = 'GRID'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class HOMO(object):
    name = 'HOMO'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('HOMOI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class KNOT(object):
    name = 'KNOT'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('KNOTI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRIM(object):
    name = 'TRIM'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('IDTRI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class XYZ(object):
    name = 'XYZ'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCNURBS'
    dtype = [('RBID', '<i8', ()), ('NPTU', '<i8', ()), ('NPTV', '<i8', ()), ('NORU', '<i8', ()), ('NORV', '<i8', ()),
             ('NSUBU', '<i8', ()), ('NSUBV', '<i8', ()), ('NGRID', '<i8', ()), ('NXYZ', '<i8', ()),
             ('NHOMO', '<i8', ()), ('NKNOT', '<i8', ()), ('NTRIM', '<i8', ()), ('POS_GRID', '<i8', ()),
             ('POS_XYZ', '<i8', ()), ('POS_HOMO', '<i8', ()), ('POS_KNOT', '<i8', ()), ('POS_TRIM', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCNURBS/GRID', '/NASTRAN/INPUT/CONTACT/BCNURBS/HOMO',
                 '/NASTRAN/INPUT/CONTACT/BCNURBS/KNOT', '/NASTRAN/INPUT/CONTACT/BCNURBS/TRIM',
                 '/NASTRAN/INPUT/CONTACT/BCNURBS/XYZ']


@register_table
class MASTERS(object):
    name = 'MASTERS'
    path = '/NASTRAN/INPUT/CONTACT/BCONECT'
    dtype = [('IDMA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SLAVES(object):
    name = 'SLAVES'
    path = '/NASTRAN/INPUT/CONTACT/BCONECT'
    dtype = [('IDSL', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCONECT'
    dtype = [('NWD', '<i8', ()), ('ID', '<i8', ()), ('BCGPID', '<i8', ()), ('BCPPID', '<i8', ()),
             ('IDSLAVE', '<i8', ()), ('IDMASTER', '<i8', ()), ('NSLA', '<i8', ()), ('NMAS', '<i8', ()),
             ('POS_SLA', '<i8', ()), ('POS_MAS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCONECT/MASTERS', '/NASTRAN/INPUT/CONTACT/BCONECT/SLAVES']


@register_table
class BCONPRG(object):
    name = 'BCONPRG'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('ID', '<i8', ()), ('U_AUGDIST', '<i8', ()), ('AUGDIST', '<f8', ()), ('U_BIAS', '<i8', ()),
             ('BIAS', '<f8', ()), ('U_CINTERF', '<i8', ()), ('CINTERF', '<f8', ()), ('U_COPTS', '<i8', ()),
             ('COPTS', '<i8', ()), ('U_COPTM', '<i8', ()), ('COPTM', '<i8', ()), ('U_ERROR', '<i8', ()),
             ('ERROR', '<f8', ()), ('U_HARDS', '<i8', ()), ('HARDS', '<f8', ()), ('U_ICOORD', '<i8', ()),
             ('ICOORD', '<i8', ()), ('U_IGLUE', '<i8', ()), ('IGLUE', '<i8', ()), ('U_ISEARCH', '<i8', ()),
             ('ISEARCH', '<i8', ()), ('U_JGLUE', '<i8', ()), ('JGLUE', '<i8', ()), ('U_PENALT', '<i8', ()),
             ('PENALT', '<f8', ()), ('U_SLIDE', '<i8', ()), ('SLIDE', '<f8', ()), ('U_STKSLP', '<i8', ()),
             ('STKSLP', '<f8', ()), ('U_TPENALT', '<i8', ()), ('TPENALT', '<f8', ()), ('U_OPINTRF', '<i8', ()),
             ('OPINTRF', '<i8', ()), ('U_TBINTRF', '<i8', ()), ('TBINTRF', '<i8', ()), ('U_CBINTRF', '<i8', ()),
             ('CBINTRF', '<i8', ()), ('U_CDINTRF', '<i8', ()), ('CDINTRF', '<i8', ()), ('U_VXINTRF', '<i8', ()),
             ('VXINTRF', '<f8', ()), ('U_VYINTRF', '<i8', ()), ('VYINTRF', '<f8', ()), ('U_VZINTRF', '<i8', ()),
             ('VZINTRF', '<f8', ()), ('U_PTINTRF', '<i8', ()), ('PTINTRF', '<f8', ()), ('U_XCINTRF', '<i8', ()),
             ('XCINTRF', '<f8', ()), ('U_YCINTRF', '<i8', ()), ('YCINTRF', '<f8', ()), ('U_ZCINTRF', '<i8', ()),
             ('ZCINTRF', '<f8', ()), ('U_OPINGP', '<i8', ()), ('OPINGP', '<i8', ()), ('U_TOLINGP', '<i8', ()),
             ('TOLINGP', '<f8', ()), ('U_CDINGP', '<i8', ()), ('CDINGP', '<i8', ()), ('U_MGINGP', '<i8', ()),
             ('MGINGP', '<f8', ()), ('U_FGCFLG', '<i8', ()), ('FGCFLG', '<i8', ()), ('U_FGCNST', '<i8', ()),
             ('FGCNST', '<f8', ()), ('U_FGCTST', '<i8', ()), ('FGCTST', '<f8', ()), ('U_FGCNSTR', '<i8', ()),
             ('FGCNSTR', '<i8', ()), ('U_FGCTSTR', '<i8', ()), ('FGCTSTR', '<i8', ()), ('U_FGCNSTI', '<i8', ()),
             ('FGCNSTI', '<i8', ()), ('U_FGCTSTI', '<i8', ()), ('FGCTSTI', '<i8', ()), ('U_FGCRCEN', '<i8', ()),
             ('FGCRCEN', '<i8', ()), ('U_SFNPNLT', '<i8', ()), ('SFNPNLT', '<f8', ()), ('U_SFTPNLT', '<i8', ()),
             ('SFTPNLT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BCONPRP(object):
    name = 'BCONPRP'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('ID', '<i8', ()), ('U_BGM', '<i8', ()), ('BGM', '<f8', ()), ('U_BGN', '<i8', ()), ('BGN', '<f8', ()),
             ('U_BGSN', '<i8', ()), ('BGSN', '<f8', ()), ('U_BGST', '<i8', ()), ('BGST', '<f8', ()),
             ('U_BNC', '<i8', ()), ('BNCI', '<i8', ()), ('BNCR', '<f8', ()), ('U_BNL', '<i8', ()), ('BNLI', '<i8', ()),
             ('BNLR', '<f8', ()), ('U_DQNEAR', '<i8', ()), ('DQNEAR', '<f8', ()), ('U_EMISS', '<i8', ()),
             ('EMISSI', '<i8', ()), ('EMISSR', '<f8', ()), ('U_FNTOL', '<i8', ()), ('FNTOL', '<f8', ()),
             ('U_FRIC', '<i8', ()), ('FRICI', '<i8', ()), ('FRICR', '<f8', ()), ('U_FRLIM', '<i8', ()),
             ('FRLIM', '<f8', ()), ('U_HBL', '<i8', ()), ('HBLI', '<i8', ()), ('HBLR', '<f8', ()), ('U_HCT', '<i8', ()),
             ('HCTI', '<i8', ()), ('HCTR', '<f8', ()), ('U_HCV', '<i8', ()), ('HCVI', '<i8', ()), ('HCVR', '<f8', ()),
             ('U_HGLUE', '<i8', ()), ('HGLUE', '<i8', ()), ('U_HNC', '<i8', ()), ('HNCI', '<i8', ()),
             ('HNCR', '<f8', ()), ('U_HNL', '<i8', ()), ('HNLI', '<i8', ()), ('HNLR', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDATA(object):
    name = 'CDATA'
    path = '/NASTRAN/INPUT/CONTACT/BCONUDS'
    dtype = [('CDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDATA(object):
    name = 'IDATA'
    path = '/NASTRAN/INPUT/CONTACT/BCONUDS'
    dtype = [('IDATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class RDATA(object):
    name = 'RDATA'
    path = '/NASTRAN/INPUT/CONTACT/BCONUDS'
    dtype = [('RDATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCONUDS'
    dtype = [('NWD', '<i8', ()), ('BID', '<i8', ()), ('BTYPE', 'S8', ()), ('GROUP', 'S8', ()), ('UNAME', 'S8', ()),
             ('NINTG', '<i8', ()), ('NREAL', '<i8', ()), ('NCHAR', '<i8', ()), ('IPOS', '<i8', ()), ('RPOS', '<i8', ()),
             ('CPOS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCONUDS/CDATA', '/NASTRAN/INPUT/CONTACT/BCONUDS/IDATA',
                 '/NASTRAN/INPUT/CONTACT/BCONUDS/RDATA']


@register_table
class BCPARA(object):
    name = 'BCPARA'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('ID', '<i8', ()), ('NBODIES', '<i8', ()), ('MAXENT', '<i8', ()), ('MAXNOD', '<i8', ()),
             ('ERROR', '<f8', ()), ('BIAS', '<f8', ()), ('ISPLIT', '<i8', ()), ('FNTOL', '<f8', ()),
             ('MAXSEP', '<i8', ()), ('ICHECK', '<i8', ()), ('ICSEP', '<i8', ()), ('IBSEP', '<i8', ()),
             ('ISHELL', '<i8', ()), ('IPRINT', '<i8', ()), ('RVCNST', '<f8', ()), ('FTYPE', '<i8', ()),
             ('FKIND', '<i8', ()), ('BEAMB', '<i8', ()), ('FSSMULT', '<f8', ()), ('FSSTOL', '<f8', ()),
             ('LINQUAD', '<i8', ()), ('INITCON', '<i8', ()), ('NGLUE', '<i8', ()), ('ITOPBM', '<i8', ()),
             ('ITOPSH', '<i8', ()), ('ITOPSD', '<i8', ()), ('NODSEP', '<i8', ()), ('METHOD', 'S4', ()),
             ('AUGMENT', '<i8', ()), ('PENALT', '<f8', ()), ('AUGDIST', '<f8', ()), ('SLDLMT', '<f8', ()),
             ('SEGSYM', '<i8', ()), ('THKOFF', '<i8', ()), ('ERRBAS', '<i8', ()), ('DDULMT', '<f8', ()),
             ('TAUGMNT', '<i8', ()), ('TPENALT', '<f8', ()), ('STKSLP', '<f8', ()), ('LINCNT', '<i8', ()),
             ('IMULTI', '<i8', ()), ('FGCFLG', '<i8', ()), ('FGCNST', '<f8', ()), ('FGCTST', '<f8', ()),
             ('FGCNSTR', '<i8', ()), ('FGCTSTR', '<i8', ()), ('FGCNSTI', '<i8', ()), ('FGCTSTI', '<i8', ()),
             ('FGCRCEN', '<i8', ()), ('SFNPNLT', '<f8', ()), ('SFTPNLT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PATCH(object):
    name = 'PATCH'
    path = '/NASTRAN/INPUT/CONTACT/BCPATCH'
    dtype = [('IDP', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()), ('G4', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCPATCH'
    dtype = [('RBID', '<i8', ()), ('NIDP', '<i8', ()), ('PATCH_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCPATCH/PATCH']


@register_table
class SEGMENT(object):
    name = 'SEGMENT'
    path = '/NASTRAN/INPUT/CONTACT/BCPFLG'
    dtype = [('IBRANCH', '<i8', ()), ('IOUTIN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCPFLG'
    dtype = [('PID', '<i8', ()), ('SEGMENT_POS', '<i8', ()), ('SEGMENT_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCPFLG/SEGMENT']


@register_table
class PIDS(object):
    name = 'PIDS'
    path = '/NASTRAN/INPUT/CONTACT/BCPROP'
    dtype = [('IP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCPROP'
    dtype = [('ID', '<i8', ()), ('POS_PIDS', '<i8', ()), ('NUM_PIDS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCPROP/PIDS']


@register_table
class RBIDS(object):
    name = 'RBIDS'
    path = '/NASTRAN/INPUT/CONTACT/BCRGSRF'
    dtype = [('RBID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCRGSRF'
    dtype = [('ID', '<i8', ()), ('NRBID', '<i8', ()), ('POS_RBID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCRGSRF/RBIDS']


@register_table
class BCRIGID(object):
    name = 'BCRIGID'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('BCRGID', '<i8', ()), ('CGID', '<i8', ()), ('CONTROL', '<i8', ()), ('NLOAD', '<i8', ()),
             ('ANGVELI', '<i8', ()), ('ANGVELR', '<f8', ()), ('DCOS1I', '<i8', ()), ('DCOS1R', '<f8', ()),
             ('DCOS2I', '<i8', ()), ('DCOS2R', '<f8', ()), ('DCOS3I', '<i8', ()), ('DCOS3R', '<f8', ()),
             ('VELRB1I', '<i8', ()), ('VELRB1R', '<f8', ()), ('VELRB2I', '<i8', ()), ('VELRB2R', '<f8', ()),
             ('VELRB3I', '<i8', ()), ('VELRB3R', '<f8', ()), ('OPTION', '<i8', ()), ('APPROV', '<i8', ()),
             ('A', '<f8', ()), ('N1', '<f8', ()), ('N2', '<f8', ()), ('N3', '<f8', ()), ('V1', '<f8', ()),
             ('V2', '<f8', ()), ('V3', '<f8', ()), ('GROW', '<i8', ()), ('GF1', '<f8', ()), ('GF2', '<f8', ()),
             ('GF3', '<f8', ()), ('TABGF1', '<i8', ()), ('TABGF2', '<i8', ()), ('TABGF3', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BCSCAP(object):
    name = 'BCSCAP'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('EID', '<i8', ()), ('IESCAP', '<i8', ()), ('NSEG', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/CONTACT/BCTABL1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCTABL1'
    dtype = [('BCID', '<i8', ()), ('NID', '<i8', ()), ('POS_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCTABL1/IDLIST']


@register_table
class MASTERS(object):
    name = 'MASTERS'
    path = '/NASTRAN/INPUT/CONTACT/BCTABLE'
    dtype = [('IDMA1', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SLAVES(object):
    name = 'SLAVES'
    path = '/NASTRAN/INPUT/CONTACT/BCTABLE'
    dtype = [('SLAVE', 'S8', ()), ('IDSLA1', '<i8', ()), ('ERROR', '<f8', ()), ('FNTOL', '<f8', ()),
             ('FRICI', '<i8', ()), ('FRICR', '<f8', ()), ('CINTERF', '<f8', ()), ('IGLUE', '<i8', ()),
             ('ISEARCH', '<i8', ()), ('ICOORD', '<i8', ()), ('JGLUE', '<i8', ()), ('TOLID', '<i8', ()),
             ('DQNEAR', '<f8', ()), ('DISTID', '<i8', ()), ('FRLIM', '<f8', ()), ('BIAS', '<f8', ()),
             ('SLIDE', '<f8', ()), ('HARDS', '<f8', ()), ('COPTS1', '<i8', ()), ('COPTM1', '<i8', ()),
             ('BKGL', 'S4', ()), ('BGSN', '<f8', ()), ('BGST', '<f8', ()), ('BGM', '<f8', ()), ('BGN', '<f8', ()),
             ('SEGS', 'S4', ()), ('PENALT', '<f8', ()), ('AUGDIST', '<f8', ()), ('TPENALT', '<f8', ()),
             ('STKSLP', '<f8', ()), ('HHHB', 'S4', ()), ('HCTI', '<i8', ()), ('HCTR', '<f8', ()), ('HCVI', '<i8', ()),
             ('HCVR', '<f8', ()), ('HNCI', '<i8', ()), ('HNCR', '<f8', ()), ('BNCI', '<i8', ()), ('BNCR', '<f8', ()),
             ('EMISSI', '<i8', ()), ('EMISSR', '<f8', ()), ('HBLI', '<i8', ()), ('HBLR', '<f8', ()),
             ('HNLI', '<i8', ()), ('HNLR', '<f8', ()), ('BNLI', '<i8', ()), ('BNLR', '<f8', ()), ('HGLUE', '<i8', ()),
             ('MASTERS', 'S8', ()), ('MASTERS_POS', '<i8', ()), ('MASTERS_LEN', '<i8', ()), ('SOL700_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SOL700(object):
    name = 'SOL700'
    path = '/NASTRAN/INPUT/CONTACT/BCTABLE'
    dtype = [('FK', '<f8', ()), ('EXP', '<f8', ()), ('METHOD', 'S4', ()), ('ADAPT', 'S4', ()), ('THICK', '<f8', ()),
             ('THICKOF', '<f8', ()), ('PENV', '<f8', ()), ('FACT', '<f8', ()), ('TSTART', '<f8', ()),
             ('TEND', '<f8', ()), ('MAXPAR', '<f8', ()), ('EROSOP', '<i8', ()), ('IADJ', '<i8', ()),
             ('SOFT', '<i8', ()), ('DEPTH', '<i8', ()), ('BSORT', '<i8', ()), ('FRCFRQ', '<i8', ()),
             ('SNLOG', '<i8', ()), ('ISYM', '<i8', ()), ('I2D3D', '<i8', ()), ('IGNORE', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCTABLE'
    dtype = [('ID', '<i8', ()), ('IDSLAVE', '<i8', ()), ('IDMAST', '<i8', ()), ('NGROUP', '<i8', ()),
             ('COPTS', '<i8', ()), ('COPTM', '<i8', ()), ('SOLNUM', '<i8', ()), ('SLAVES_POS', '<i8', ()),
             ('SLAVES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCTABLE/MASTERS', '/NASTRAN/INPUT/CONTACT/BCTABLE/SLAVES',
                 '/NASTRAN/INPUT/CONTACT/BCTABLE/SOL700']


@register_table
class HOMO(object):
    name = 'HOMO'
    path = '/NASTRAN/INPUT/CONTACT/BCTRIM'
    dtype = [('HOMOI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class ISOP(object):
    name = 'ISOP'
    path = '/NASTRAN/INPUT/CONTACT/BCTRIM'
    dtype = [('XISOPARM', '<f8', ()), ('YISOPARM', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class KNOT(object):
    name = 'KNOT'
    path = '/NASTRAN/INPUT/CONTACT/BCTRIM'
    dtype = [('KNOTI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BCTRIM'
    dtype = [('ID', '<i8', ()), ('NPTU', '<i8', ()), ('NORU', '<i8', ()), ('NSUB', '<i8', ()), ('NISOP', '<i8', ()),
             ('NHOMO', '<i8', ()), ('NKNOT', '<i8', ()), ('POS_ISO', '<i8', ()), ('POS_HOMO', '<i8', ()),
             ('POS_KNOT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BCTRIM/HOMO', '/NASTRAN/INPUT/CONTACT/BCTRIM/ISOP',
                 '/NASTRAN/INPUT/CONTACT/BCTRIM/KNOT']


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/CONTACT/BLSEG'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BLSEG'
    dtype = [('ID', '<i8', ()), ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BLSEG/GRIDS']


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/CONTACT/BOUTPUT'
    dtype = [('G1', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BOUTPUT'
    dtype = [('ID', '<i8', ()), ('POS_GRIDS', '<i8', ()), ('NUM_GRIDS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BOUTPUT/GRIDS']


@register_table
class BSQUEAL(object):
    name = 'BSQUEAL'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('ID', '<i8', ()), ('OMETH', '<f8', ()), ('AVSTIF', '<f8', ()), ('BSONLY', 'S4', ()), ('R1', '<f8', ()),
             ('R2', '<f8', ()), ('R3', '<f8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/CONTACT/BSURF'
    dtype = [('IDI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRUBY(object):
    name = 'THRUBY'
    path = '/NASTRAN/INPUT/CONTACT/BSURF'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('INC', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BSURF'
    dtype = [('ID', '<i8', ()), ('DOMAIN_ID', '<i8', ()), ('NUM_IDS', '<i8', ()), ('NUM_THRUBY', '<i8', ()),
             ('POS_IDS', '<i8', ()), ('POS_THRUBY', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BSURF/IDS', '/NASTRAN/INPUT/CONTACT/BSURF/THRUBY']


@register_table
class LIST(object):
    name = 'LIST'
    path = '/NASTRAN/INPUT/CONTACT/BSURF_OLD'
    dtype = [('EID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THBY(object):
    name = 'THBY'
    path = '/NASTRAN/INPUT/CONTACT/BSURF_OLD'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('INC', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/BSURF_OLD'
    dtype = [('ID', '<i8', ()), ('LIST_POS', '<i8', ()), ('LIST_LEN', '<i8', ()), ('THBY_POS', '<i8', ()),
             ('THBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/BSURF_OLD/LIST', '/NASTRAN/INPUT/CONTACT/BSURF_OLD/THBY']


@register_table
class PRJCON(object):
    name = 'PRJCON'
    path = '/NASTRAN/INPUT/CONTACT'
    dtype = [('BID', '<i8', ()), ('HEAT', '<i8', ()), ('SET3M1', '<i8', ()), ('SET3S1', '<i8', ()), ('H', '<f8', ()),
             ('PCONID', '<i8', ()), ('F3', '<f8', ()), ('EMISM', '<f8', ()), ('EMISS', '<f8', ()), ('F4', '<f8', ()),
             ('RADCM', '<i8', ()), ('RADCS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/CONTACT/UNGLUE'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/CONTACT/UNGLUE'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRUBY(object):
    name = 'THRUBY'
    path = '/NASTRAN/INPUT/CONTACT/UNGLUE'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('STEP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/CONTACT/UNGLUE'
    dtype = [('ID', '<i8', ()), ('BID', '<i8', ()), ('IDLIST_POS', '<i8', ()), ('IDLIST_LEN', '<i8', ()),
             ('THRU_POS', '<i8', ()), ('THRU_LEN', '<i8', ()), ('THRUBY_POS', '<i8', ()), ('THRUBY_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/CONTACT/UNGLUE/IDLIST', '/NASTRAN/INPUT/CONTACT/UNGLUE/THRU',
                 '/NASTRAN/INPUT/CONTACT/UNGLUE/THRUBY']


@register_table
class CORD1C(object):
    name = 'CORD1C'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD1R(object):
    name = 'CORD1R'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD1S(object):
    name = 'CORD1S'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD2C(object):
    name = 'CORD2C'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('RID', '<i8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()), ('A3', '<f8', ()),
             ('B1', '<f8', ()), ('B2', '<f8', ()), ('B3', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()),
             ('C3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD2R(object):
    name = 'CORD2R'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('RID', '<i8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()), ('A3', '<f8', ()),
             ('B1', '<f8', ()), ('B2', '<f8', ()), ('B3', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()),
             ('C3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD2S(object):
    name = 'CORD2S'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('RID', '<i8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()), ('A3', '<f8', ()),
             ('B1', '<f8', ()), ('B2', '<f8', ()), ('B3', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()),
             ('C3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD3G(object):
    name = 'CORD3G'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('METHOD', 'S8', ()), ('FORM', 'S8', ()), ('THETAID', '<i8', (3,)),
             ('CIDREF', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORD3R(object):
    name = 'CORD3R'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM'
    dtype = [('CID', '<i8', ()), ('N1', '<i8', ()), ('N2', '<i8', ()), ('N3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDATA(object):
    name = 'IDATA'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION'
    dtype = [('DATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class RDATA(object):
    name = 'RDATA'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION'
    dtype = [('DATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION'
    dtype = [('CID', '<i8', ()), ('TYPE', '<i8', ()), ('IINDEX', '<i8', ()), ('RINDEX', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION/IDATA',
                 '/NASTRAN/INPUT/COORDINATE_SYSTEM/TRANSFORMATION/RDATA']


@register_table
class DIRECTION(object):
    name = 'DIRECTION'
    path = '/NASTRAN/INPUT/DESIGN/BEADVAR'
    dtype = [('NORM', '<i8', ()), ('CID1', '<i8', ()), ('N', '<f8', (3,)), ('XLB', '<f8', ()), ('XUB', '<f8', ()),
             ('DELXV', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GRIDSET(object):
    name = 'GRIDSET'
    path = '/NASTRAN/INPUT/DESIGN/BEADVAR'
    dtype = [('NDGSET', '<i8', ()), ('DGSET', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PATTERN(object):
    name = 'PATTERN'
    path = '/NASTRAN/INPUT/DESIGN/BEADVAR'
    dtype = [('CID3', '<i8', ()), ('PATTY', '<i8', ()), ('PATDIR', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SYMPLANE(object):
    name = 'SYMPLANE'
    path = '/NASTRAN/INPUT/DESIGN/BEADVAR'
    dtype = [('CID2', '<i8', ()), ('MS1', '<i8', ()), ('MS2', '<i8', ()), ('MS3', '<i8', ()), ('CSM', '<i8', ()),
             ('NCSM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/BEADVAR'
    dtype = [('BEADID', '<i8', ()), ('PTYPE', 'S8', ()), ('PID', '<i8', ()), ('MW', '<f8', ()), ('MH', '<f8', ()),
             ('ANG', '<f8', ()), ('BF', '<i8', ()), ('SKIP', '<i8', ()), ('DIR_LEN', '<i8', ()), ('SYM_LEN', '<i8', ()),
             ('PAT_LEN', '<i8', ()), ('GRID_LEN', '<i8', ()), ('DIR_POS', '<i8', ()), ('SYM_POS', '<i8', ()),
             ('PAT_POS', '<i8', ()), ('GRID_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/BEADVAR/DIRECTION', '/NASTRAN/INPUT/DESIGN/BEADVAR/GRIDSET',
                 '/NASTRAN/INPUT/DESIGN/BEADVAR/PATTERN', '/NASTRAN/INPUT/DESIGN/BEADVAR/SYMPLANE']


@register_table
class DCLIST(object):
    name = 'DCLIST'
    path = '/NASTRAN/INPUT/DESIGN/DCONADD'
    dtype = [('DC', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DCONADD'
    dtype = [('DCID', '<i8', ()), ('DCLIST_POS', '<i8', ()), ('DCLIST_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DCONADD/DCLIST']


@register_table
class DCONSTR(object):
    name = 'DCONSTR'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('DCID', '<i8', ()), ('RID', '<i8', ()), ('LALLOW', '<f8', ()), ('UALLOW', '<f8', ()),
             ('LOWFQ', '<f8', ()), ('HIGHFQ', '<f8', ()), ('DTYPE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DVAL(object):
    name = 'DVAL'
    path = '/NASTRAN/INPUT/DESIGN/DDVAL'
    dtype = [('V', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DDVAL'
    dtype = [('DVID', '<i8', ()), ('DVAL_POS', '<i8', ()), ('DVAL_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DDVAL/DVAL']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('XINIT', '<f8', ()), ('XLB', '<f8', ()), ('XUB', '<f8', ()),
             ('DELX', '<f8', ()), ('DVID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDVC(object):
    name = 'IDVC'
    path = '/NASTRAN/INPUT/DESIGN/DLINK'
    dtype = [('INDV', '<i8', ()), ('C', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DLINK'
    dtype = [('ID', '<i8', ()), ('DVID', '<i8', ()), ('C0', '<f8', ()), ('CMULT', '<f8', ()), ('IDVC_POS', '<i8', ()),
             ('IDVC_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DLINK/IDVC']


@register_table
class DOPTPRM(object):
    name = 'DOPTPRM'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('APRCOD', '<i8', ()), ('IPRINT', '<i8', ()), ('DESMAX', '<i8', ()), ('METHOD', '<i8', ()),
             ('DELP', '<f8', ()), ('DPMIN', '<f8', ()), ('PTOL', '<f8', ()), ('CONV1', '<f8', ()), ('CONV2', '<f8', ()),
             ('GMAX', '<f8', ()), ('DELX', '<f8', ()), ('DXMIN', '<f8', ()), ('DELB', '<f8', ()), ('GSCAL', '<f8', ()),
             ('CONVDV', '<f8', ()), ('CONVPR', '<f8', ()), ('P1', '<i8', ()), ('P2', '<i8', ()), ('CT', '<f8', ()),
             ('CTMIN', '<f8', ()), ('DABOBJ', '<f8', ()), ('DELOBJ', '<f8', ()), ('DOBJ1', '<f8', ()),
             ('DOBJ2', '<f8', ()), ('DX1', '<f8', ()), ('DX2', '<f8', ()), ('ISCAL', '<i8', ()), ('ITMAX', '<i8', ()),
             ('ITRMOP', '<i8', ()), ('IWRITE', '<i8', ()), ('IGMAX', '<i8', ()), ('JTMAX', '<i8', ()),
             ('ITRMST', '<i8', ()), ('JPRINT', '<i8', ()), ('IPRNT1', '<i8', ()), ('IPRNT2', '<i8', ()),
             ('JWRITE', '<i8', ()), ('STPSCL', '<f8', ()), ('FSDMAX', '<i8', ()), ('FSDALP', '<f8', ()),
             ('DISCOD', '<i8', ()), ('DISBEG', '<i8', ()), ('PLVIOL', '<i8', ()), ('P2CR', '<i8', ()),
             ('P2CDDV', '<i8', ()), ('P2CP', '<i8', ()), ('P2CC', '<i8', ()), ('P2CM', '<i8', ()), ('P2CBL', '<i8', ()),
             ('P2RSET', '<i8', ()), ('P2CALL', '<i8', ()), ('ADSCOD', '<i8', ()), ('PENAL', '<f8', ()),
             ('DRATIO', '<f8', ()), ('AUTOSE', '<i8', ()), ('NUMDIV', '<i8', ()), ('TCHECK', '<i8', ()),
             ('TDMIN', '<f8', ()), ('TFWET', '<f8', ()), ('TREGION', '<i8', ()), ('EATA1', '<f8', ()),
             ('EATA2', '<f8', ()), ('EATA3', '<f8', ()), ('UPDFAC1', '<f8', ()), ('UPDFAC2', '<f8', ()),
             ('DPMAX', '<f8', ()), ('DXMAX', '<f8', ()), ('OPTCOD', 'S8', ()), ('OBJMOD', '<i8', ()),
             ('DELXESL', '<f8', ()), ('DSMXESL', '<i8', ()), ('NASPR0', '<i8', ()), ('LSCOD', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ATTI(object):
    name = 'ATTI'
    path = '/NASTRAN/INPUT/DESIGN/DRESP1'
    dtype = [('ATTI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DRESP1'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('RTYPE', 'S8', ()), ('PTYPE', 'S8', ()), ('RPSID', '<i8', ()),
             ('NTUSED', 'S4', ()), ('AFPMID', '<i8', ()), ('REGION', '<i8', ()), ('ATTA', '<i8', ()),
             ('ATTBI', '<i8', ()), ('ATTBR', '<f8', ()), ('ATTI_LEN', '<i8', ()), ('ATTI_POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DRESP1/ATTI']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DNODE(object):
    name = 'DNODE'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('GRID', '<i8', ()), ('COMP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP1(object):
    name = 'DRESP1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('R1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP2(object):
    name = 'DRESP2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('R2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('LABEL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL1(object):
    name = 'DVCREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DC1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL2(object):
    name = 'DVCREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DC2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL1(object):
    name = 'DVMREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DM1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL2(object):
    name = 'DVMREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DM2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL1(object):
    name = 'DVPREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DP1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL2(object):
    name = 'DVPREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('DP2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DRESP2'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('EQID', '<i8', ()), ('REGION', '<i8', ()), ('METHOD', '<i8', ()),
             ('C1', '<f8', ()), ('C2', '<f8', ()), ('C3', '<f8', ()), ('DESVAR_POS', '<i8', ()),
             ('DTABLE_POS', '<i8', ()), ('DRESP1_POS', '<i8', ()), ('DNODE_POS', '<i8', ()), ('DVPREL1_POS', '<i8', ()),
             ('DVCREL1_POS', '<i8', ()), ('DVMREL1_POS', '<i8', ()), ('DVPREL2_POS', '<i8', ()),
             ('DVCREL2_POS', '<i8', ()), ('DVMREL2_POS', '<i8', ()), ('DRESP2_POS', '<i8', ()),
             ('DESVAR_LEN', '<i8', ()), ('DTABLE_LEN', '<i8', ()), ('DRESP1_LEN', '<i8', ()), ('DNODE_LEN', '<i8', ()),
             ('DVPREL1_LEN', '<i8', ()), ('DVCREL1_LEN', '<i8', ()), ('DVMREL1_LEN', '<i8', ()),
             ('DVPREL2_LEN', '<i8', ()), ('DVCREL2_LEN', '<i8', ()), ('DVMREL2_LEN', '<i8', ()),
             ('DRESP2_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DRESP2/DESVAR', '/NASTRAN/INPUT/DESIGN/DRESP2/DNODE',
                 '/NASTRAN/INPUT/DESIGN/DRESP2/DRESP1', '/NASTRAN/INPUT/DESIGN/DRESP2/DRESP2',
                 '/NASTRAN/INPUT/DESIGN/DRESP2/DTABLE', '/NASTRAN/INPUT/DESIGN/DRESP2/DVCREL1',
                 '/NASTRAN/INPUT/DESIGN/DRESP2/DVCREL2', '/NASTRAN/INPUT/DESIGN/DRESP2/DVMREL1',
                 '/NASTRAN/INPUT/DESIGN/DRESP2/DVMREL2', '/NASTRAN/INPUT/DESIGN/DRESP2/DVPREL1',
                 '/NASTRAN/INPUT/DESIGN/DRESP2/DVPREL2']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DNODE(object):
    name = 'DNODE'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('GRID', '<i8', ()), ('COMP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP1(object):
    name = 'DRESP1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('R1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP2(object):
    name = 'DRESP2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('R2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('LABE', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL1(object):
    name = 'DVCREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DC1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL2(object):
    name = 'DVCREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DC2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVLREL1(object):
    name = 'DVLREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DL1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL1(object):
    name = 'DVMREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DM1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL2(object):
    name = 'DVMREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DM2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL1(object):
    name = 'DVPREL1'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DP1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL2(object):
    name = 'DVPREL2'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('DP2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GRDTYP(object):
    name = 'GRDTYP'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('GRDTYP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class USRDATA(object):
    name = 'USRDATA'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('USRDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DRESP3'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('GRPNM', 'S8', ()), ('TYPNM', 'S8', ()), ('REGION', '<i8', ()),
             ('NR3', '<i8', ()), ('DESVAR_POS', '<i8', ()), ('DTABLE_POS', '<i8', ()), ('DRESP1_POS', '<i8', ()),
             ('DNODE_POS', '<i8', ()), ('DVPREL1_POS', '<i8', ()), ('DVCREL1_POS', '<i8', ()),
             ('DVMREL1_POS', '<i8', ()), ('DVPREL2_POS', '<i8', ()), ('DVCREL2_POS', '<i8', ()),
             ('DVMREL2_POS', '<i8', ()), ('DRESP2_POS', '<i8', ()), ('DVLREL1_POS', '<i8', ()),
             ('USRDATA_POS', '<i8', ()), ('GRDTYP_POS', '<i8', ()), ('DESVAR_LEN', '<i8', ()),
             ('DTABLE_LEN', '<i8', ()), ('DRESP1_LEN', '<i8', ()), ('DNODE_LEN', '<i8', ()), ('DVPREL1_LEN', '<i8', ()),
             ('DVCREL1_LEN', '<i8', ()), ('DVMREL1_LEN', '<i8', ()), ('DVPREL2_LEN', '<i8', ()),
             ('DVCREL2_LEN', '<i8', ()), ('DVMREL2_LEN', '<i8', ()), ('DRESP2_LEN', '<i8', ()),
             ('DVLREL1_LEN', '<i8', ()), ('USRDATA_LEN', '<i8', ()), ('GRDTYP_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DRESP3/DESVAR', '/NASTRAN/INPUT/DESIGN/DRESP3/DNODE',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/DRESP1', '/NASTRAN/INPUT/DESIGN/DRESP3/DRESP2',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/DTABLE', '/NASTRAN/INPUT/DESIGN/DRESP3/DVCREL1',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/DVCREL2', '/NASTRAN/INPUT/DESIGN/DRESP3/DVLREL1',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/DVMREL1', '/NASTRAN/INPUT/DESIGN/DRESP3/DVMREL2',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/DVPREL1', '/NASTRAN/INPUT/DESIGN/DRESP3/DVPREL2',
                 '/NASTRAN/INPUT/DESIGN/DRESP3/GRDTYP', '/NASTRAN/INPUT/DESIGN/DRESP3/USRDATA']


@register_table
class DSCREEN(object):
    name = 'DSCREEN'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('TYPE', '<i8', ()), ('TRS', '<f8', ()), ('NSTR', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VALUES(object):
    name = 'VALUES'
    path = '/NASTRAN/INPUT/DESIGN/DTABLE'
    dtype = [('LABEL', 'S8', ()), ('VALUE', '<f8', ()), ('VALUEI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DTABLE'
    dtype = [('DOMAIN_ID', '<i8', ()), ('START', '<i8', ()), ('LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DTABLE/VALUES']


@register_table
class VALUES(object):
    name = 'VALUES'
    path = '/NASTRAN/INPUT/DESIGN/DTABLE2'
    dtype = [('LABEL', 'S8', ()), ('PNAME', 'S8', ()), ('PID', '<i8', ()), ('FNAME', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DTABLE2'
    dtype = [('DOMAIN_ID', '<i8', ()), ('START', '<i8', ()), ('LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DTABLE2/VALUES']


@register_table
class TERMS(object):
    name = 'TERMS'
    path = '/NASTRAN/INPUT/DESIGN/DVCREL1'
    dtype = [('DVID', '<i8', ()), ('COEF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVCREL1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('EID', '<i8', ()), ('FID', '<i8', ()), ('CPMIN', '<f8', ()),
             ('CPMAX', '<f8', ()), ('C0', '<f8', ()), ('CPNAME', 'S8', ()), ('TERMS_POS', '<i8', ()),
             ('TERMS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVCREL1/TERMS']


@register_table
class DVIDS(object):
    name = 'DVIDS'
    path = '/NASTRAN/INPUT/DESIGN/DVCREL2'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class LABLS(object):
    name = 'LABLS'
    path = '/NASTRAN/INPUT/DESIGN/DVCREL2'
    dtype = [('LABL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVCREL2'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('EID', '<i8', ()), ('FID', '<i8', ()), ('CPMIN', '<f8', ()),
             ('CPMAX', '<f8', ()), ('EQID', '<i8', ()), ('CPNAME', 'S8', ()), ('DVIDS_POS', '<i8', ()),
             ('DVIDS_LEN', '<i8', ()), ('LABLS_POS', '<i8', ()), ('LABLS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVCREL2/DVIDS', '/NASTRAN/INPUT/DESIGN/DVCREL2/LABLS']


@register_table
class DVGRID(object):
    name = 'DVGRID'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('DVID', '<i8', ()), ('NODE', '<i8', ()), ('CID', '<i8', ()), ('COEFF', '<f8', ()), ('N', '<f8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COEFS(object):
    name = 'COEFS'
    path = '/NASTRAN/INPUT/DESIGN/DVLREL1'
    dtype = [('DVID', '<i8', ()), ('COEF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVLREL1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('SID', '<i8', ()), ('LNAME', 'S8', ()), ('LMIN', '<f8', ()),
             ('LMAX', '<f8', ()), ('C0', '<f8', ()), ('ATT1', '<i8', ()), ('ATT2', '<i8', ()), ('ATT3', '<i8', ()),
             ('ATT4', '<i8', ()), ('ATT5', '<i8', ()), ('COEFS_POS', '<i8', ()), ('COEFS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVLREL1/COEFS']


@register_table
class DVIDS(object):
    name = 'DVIDS'
    path = '/NASTRAN/INPUT/DESIGN/DVLREL2'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class LABLS(object):
    name = 'LABLS'
    path = '/NASTRAN/INPUT/DESIGN/DVLREL2'
    dtype = [('LABL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVLREL2'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('SID', '<i8', ()), ('LNAME', 'S8', ()), ('LMIN', '<f8', ()),
             ('LMAX', '<f8', ()), ('EQID', '<i8', ()), ('ATT1', '<i8', ()), ('ATT2', '<i8', ()), ('ATT3', '<i8', ()),
             ('ATT4', '<i8', ()), ('ATT5', '<i8', ()), ('DVIDS_POS', '<i8', ()), ('DVIDS_LEN', '<i8', ()),
             ('LABLS_POS', '<i8', ()), ('LABLS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVLREL2/DVIDS', '/NASTRAN/INPUT/DESIGN/DVLREL2/LABLS']


@register_table
class RELATION(object):
    name = 'RELATION'
    path = '/NASTRAN/INPUT/DESIGN/DVMREL1'
    dtype = [('DVID', '<i8', ()), ('COEF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVMREL1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('MID', '<i8', ()), ('FID', '<i8', ()), ('PMIN', '<f8', ()),
             ('PMAX', '<f8', ()), ('C0', '<f8', ()), ('MPNAME', 'S8', ()), ('START', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVMREL1/RELATION']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/DVMREL2'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/DVMREL2'
    dtype = [('CID', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVMREL2'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('MID', '<i8', ()), ('FID', '<i8', ()), ('PMIN', '<f8', ()),
             ('PMAX', '<f8', ()), ('EQID', '<i8', ()), ('MPNAME', 'S8', ()), ('DESVAR_POS', '<i8', ()),
             ('DTABLE_POS', '<i8', ()), ('DESVAR_LEN', '<i8', ()), ('DTABLE_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVMREL2/DESVAR', '/NASTRAN/INPUT/DESIGN/DVMREL2/DTABLE']


@register_table
class RELATION(object):
    name = 'RELATION'
    path = '/NASTRAN/INPUT/DESIGN/DVPREL1'
    dtype = [('DVID', '<i8', ()), ('COEF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVPREL1'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('PID', '<i8', ()), ('FID', '<i8', ()), ('PMIN', '<f8', ()),
             ('PMAX', '<f8', ()), ('C0', '<f8', ()), ('PNAME', 'S8', ()), ('START', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVPREL1/RELATION']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/DVPREL2'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/DVPREL2'
    dtype = [('CID', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVPREL2'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S8', ()), ('PID', '<i8', ()), ('FID', '<i8', ()), ('PMIN', '<f8', ()),
             ('PMAX', '<f8', ()), ('EQID', '<i8', ()), ('PNAME', 'S8', ()), ('DESVAR_POS', '<i8', ()),
             ('DTABLE_POS', '<i8', ()), ('DESVAR_LEN', '<i8', ()), ('DTABLE_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVPREL2/DESVAR', '/NASTRAN/INPUT/DESIGN/DVPREL2/DTABLE']


@register_table
class DVPSURF(object):
    name = 'DVPSURF'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('SUBID', '<i8', ()), ('DVID', '<i8', ()), ('COEF', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SFACTS(object):
    name = 'SFACTS'
    path = '/NASTRAN/INPUT/DESIGN/DVSHAP'
    dtype = [('DVID', '<i8', ()), ('COLID', '<i8', ()), ('SF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/DVSHAP'
    dtype = [('SFACTS_POS', '<i8', ()), ('SFACTS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/DVSHAP/SFACTS']


@register_table
class MODTRAK(object):
    name = 'MODTRAK'
    path = '/NASTRAN/INPUT/DESIGN'
    dtype = [('SID', '<i8', ()), ('RANGELO', '<i8', ()), ('RANGEHI', '<i8', ()), ('MFILTER', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INDEP(object):
    name = 'INDEP'
    path = '/NASTRAN/INPUT/DESIGN/SEDLINK'
    dtype = [('SEID', '<i8', ()), ('INDV', '<i8', ()), ('C', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/SEDLINK'
    dtype = [('ID', '<i8', ()), ('DSEID', '<i8', ()), ('DDVID', '<i8', ()), ('C0', '<f8', ()), ('CMULT', '<f8', ()),
             ('INDEP_POS', '<i8', ()), ('INDEP_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/SEDLINK/INDEP']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DNODE(object):
    name = 'DNODE'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('GRID', '<i8', ()), ('COMP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP1(object):
    name = 'DRESP1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('R1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP2(object):
    name = 'DRESP2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('R2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('LABL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL1(object):
    name = 'DVCREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DC1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL2(object):
    name = 'DVCREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DC2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL1(object):
    name = 'DVMREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DM1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL2(object):
    name = 'DVMREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DM2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL1(object):
    name = 'DVPREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DP1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL2(object):
    name = 'DVPREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('SEID', '<i8', ()), ('DP2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP2'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('EQID', '<i8', ()), ('REGION', '<i8', ()), ('METHOD', '<i8', ()),
             ('C1', '<f8', ()), ('C2', '<f8', ()), ('C3', '<f8', ()), ('DESVAR_POS', '<i8', ()),
             ('DESVAR_LEN', '<i8', ()), ('DTABLE_POS', '<i8', ()), ('DTABLE_LEN', '<i8', ()), ('DRESP1_POS', '<i8', ()),
             ('DRESP1_LEN', '<i8', ()), ('DNODE_POS', '<i8', ()), ('DNODE_LEN', '<i8', ()), ('DVPREL1_POS', '<i8', ()),
             ('DVPREL1_LEN', '<i8', ()), ('DVCREL1_POS', '<i8', ()), ('DVCREL1_LEN', '<i8', ()),
             ('DVMREL1_POS', '<i8', ()), ('DVMREL1_LEN', '<i8', ()), ('DVPREL2_POS', '<i8', ()),
             ('DVPREL2_LEN', '<i8', ()), ('DVCREL2_POS', '<i8', ()), ('DVCREL2_LEN', '<i8', ()),
             ('DVMREL2_POS', '<i8', ()), ('DVMREL2_LEN', '<i8', ()), ('DRESP2_POS', '<i8', ()),
             ('DRESP2_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/SEDRSP2/DESVAR', '/NASTRAN/INPUT/DESIGN/SEDRSP2/DNODE',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP2/DRESP1', '/NASTRAN/INPUT/DESIGN/SEDRSP2/DRESP2',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP2/DTABLE', '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVCREL1',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVCREL2', '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVMREL1',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVMREL2', '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVPREL1',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP2/DVPREL2']


@register_table
class DESVAR(object):
    name = 'DESVAR'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DNODE(object):
    name = 'DNODE'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('GRID', '<i8', ()), ('COMP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP1(object):
    name = 'DRESP1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('R1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DRESP2(object):
    name = 'DRESP2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('R2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DTABLE(object):
    name = 'DTABLE'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('LABL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL1(object):
    name = 'DVCREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DC1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVCREL2(object):
    name = 'DVCREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DC2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVLREL1(object):
    name = 'DVLREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('L1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL1(object):
    name = 'DVMREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DM1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVMREL2(object):
    name = 'DVMREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DM2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL1(object):
    name = 'DVPREL1'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DP1ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DVPREL2(object):
    name = 'DVPREL2'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('SEID', '<i8', ()), ('DP2ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GRDTYPS(object):
    name = 'GRDTYPS'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('GRDTYP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class USER(object):
    name = 'USER'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('USRDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/SEDRSP3'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('GRPNM', 'S8', ()), ('TYPNM', 'S8', ()), ('REGION', '<i8', ()),
             ('NR3', '<i8', ()), ('DESVAR_POS', '<i8', ()), ('DESVAR_LEN', '<i8', ()), ('DTABLE_POS', '<i8', ()),
             ('DTABLE_LEN', '<i8', ()), ('DRESP1_POS', '<i8', ()), ('DRESP1_LEN', '<i8', ()), ('DNODE_POS', '<i8', ()),
             ('DNODE_LEN', '<i8', ()), ('DVPREL1_POS', '<i8', ()), ('DVPREL1_LEN', '<i8', ()),
             ('DVCREL1_POS', '<i8', ()), ('DVCREL1_LEN', '<i8', ()), ('DVMREL1_POS', '<i8', ()),
             ('DVMREL1_LEN', '<i8', ()), ('DVPREL2_POS', '<i8', ()), ('DVPREL2_LEN', '<i8', ()),
             ('DVCREL2_POS', '<i8', ()), ('DVCREL2_LEN', '<i8', ()), ('DVMREL2_POS', '<i8', ()),
             ('DVMREL2_LEN', '<i8', ()), ('DRESP2_POS', '<i8', ()), ('DRESP2_LEN', '<i8', ()),
             ('DVLREL1_POS', '<i8', ()), ('DVLREL1_LEN', '<i8', ()), ('USER_POS', '<i8', ()), ('USER_LEN', '<i8', ()),
             ('GRDTYPS_POS', '<i8', ()), ('GRDTYPS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/SEDRSP3/DESVAR', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DNODE',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/DRESP1', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DRESP2',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/DTABLE', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVCREL1',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVCREL2', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVLREL1',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVMREL1', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVMREL2',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVPREL1', '/NASTRAN/INPUT/DESIGN/SEDRSP3/DVPREL2',
                 '/NASTRAN/INPUT/DESIGN/SEDRSP3/GRDTYPS', '/NASTRAN/INPUT/DESIGN/SEDRSP3/USER']


@register_table
class DISCRETE(object):
    name = 'DISCRETE'
    path = '/NASTRAN/INPUT/DESIGN/TOMVAR'
    dtype = [('DVID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SYMPLANE(object):
    name = 'SYMPLANE'
    path = '/NASTRAN/INPUT/DESIGN/TOMVAR'
    dtype = [('CID', '<i8', ()), ('MS1', '<i8', ()), ('MS2', '<i8', ()), ('MS3', '<i8', ()), ('CSM', '<i8', ()),
             ('NCSM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TOMVARID(object):
    name = 'TOMVARID'
    path = '/NASTRAN/INPUT/DESIGN/TOMVAR'
    dtype = [('TID', '<i8', ()), ('C0', '<f8', ()), ('C1', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/TOMVAR'
    dtype = [('TOMID', '<i8', ()), ('PTYPE', 'S8', ()), ('PID', '<i8', ()), ('FID', '<i8', ()), ('PNAME', 'S8', ()),
             ('XINIT', '<f8', ()), ('XLB', '<f8', ()), ('XUB', '<f8', ()), ('DELXV', '<f8', ()), ('SYM_LEN', '<i8', ()),
             ('DIS_LEN', '<i8', ()), ('TOM_LEN', '<i8', ()), ('SYM_POS', '<i8', ()), ('DIS_POS', '<i8', ()),
             ('TOM_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/TOMVAR/DISCRETE', '/NASTRAN/INPUT/DESIGN/TOMVAR/SYMPLANE',
                 '/NASTRAN/INPUT/DESIGN/TOMVAR/TOMVARID']


@register_table
class CAST(object):
    name = 'CAST'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('CID2', '<i8', ()), ('DDI', '<i8', ()), ('DIE', '<i8', ()), ('ALIGN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class EXTRUSION(object):
    name = 'EXTRUSION'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('CID3', '<i8', ()), ('EDI', '<i8', ()), ('ALIGN2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class MIN(object):
    name = 'MIN'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('TDMIN', '<f8', ()), ('TDMAX', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class STRESS(object):
    name = 'STRESS'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('STLIM', '<f8', ()), ('PNORM', '<f8', ()), ('PSFT', '<f8', ()), ('PMETHOD', '<i8', ()),
             ('PMATCH', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SYMPLANE(object):
    name = 'SYMPLANE'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('CID', '<i8', ()), ('MS1', '<i8', ()), ('MS2', '<i8', ()), ('MS3', '<i8', ()), ('CSM', '<i8', ()),
             ('NCSM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DESIGN/TOPVAR'
    dtype = [('TOPID', '<i8', ()), ('LABEL', 'S8', ()), ('PTYPE', 'S8', ()), ('XINIT', '<f8', ()), ('XLB', '<f8', ()),
             ('DELXV', '<f8', ()), ('POWER', '<f8', ()), ('PID', '<i8', ()), ('SYM_LEN', '<i8', ()),
             ('CAST_LEN', '<i8', ()), ('EXT_LEN', '<i8', ()), ('MIN_LEN', '<i8', ()), ('STRESS_LEN', '<i8', ()),
             ('SYM_POS', '<i8', ()), ('CAST_POS', '<i8', ()), ('EXT_POS', '<i8', ()), ('MIN_POS', '<i8', ()),
             ('STRESS_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DESIGN/TOPVAR/CAST', '/NASTRAN/INPUT/DESIGN/TOPVAR/EXTRUSION',
                 '/NASTRAN/INPUT/DESIGN/TOPVAR/MIN', '/NASTRAN/INPUT/DESIGN/TOPVAR/STRESS',
                 '/NASTRAN/INPUT/DESIGN/TOPVAR/SYMPLANE']


@register_table
class DOMAINS(object):
    name = 'DOMAINS'
    path = '/NASTRAN/INPUT'
    dtype = [('ID', '<i8', ()), ('SE', '<i8', ()), ('AFPM', '<i8', ()), ('TRMC', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACLOAD(object):
    name = 'ACLOAD'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('MUNIT', '<i8', ()), ('PUNIT', '<i8', ()), ('SCLRE', '<f8', ()), ('SCLIM', '<f8', ()),
             ('LSQN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACTRIM(object):
    name = 'ACTRIM'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('NAME', 'S8', ()), ('MUNIT', '<i8', ()), ('PUNIT', '<i8', ()), ('SCLRE', '<f8', ()), ('SCLIM', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAMPBLL(object):
    name = 'CAMPBLL'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('CID', '<i8', ()), ('VPARM', '<i8', ()), ('DDVALID', '<i8', ()), ('TYPE', 'S8', ()), ('ID', '<i8', ()),
             ('FID', '<i8', ()), ('NAME', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGB(object):
    name = 'EIGB'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('METHOD', 'S8', ()), ('L1', '<f8', ()), ('L2', '<f8', ()), ('NEP', '<i8', ()),
             ('NDP', '<i8', ()), ('NDN', '<i8', ()), ('NORM', 'S8', ()), ('G', '<i8', ()), ('C', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class REGION(object):
    name = 'REGION'
    path = '/NASTRAN/INPUT/DYNAMIC/EIGC'
    dtype = [('AAJ', '<f8', ()), ('WAJ', '<f8', ()), ('ABJ', '<f8', ()), ('WBJ', '<f8', ()), ('LJ', '<f8', ()),
             ('NEJ', '<i8', ()), ('NDJ', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/EIGC'
    dtype = [('SID', '<i8', ()), ('METHOD', 'S8', ()), ('NORM', 'S8', ()), ('G', '<i8', ()), ('C', '<i8', ()),
             ('E', '<f8', ()), ('ND1', '<i8', ()), ('REGION_POS', '<i8', ()), ('REGION_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/EIGC/REGION']


@register_table
class EIGP(object):
    name = 'EIGP'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('METHOD', 'S8', ()), ('ALPHA', '<f8', ()), ('OMEGA', '<f8', ()), ('M', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGR(object):
    name = 'EIGR'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('METHOD', 'S8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('NE', '<i8', ()),
             ('ND', '<i8', ()), ('NORM', 'S8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FREQS(object):
    name = 'FREQS'
    path = '/NASTRAN/INPUT/DYNAMIC/EIGRL'
    dtype = [('FI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/EIGRL'
    dtype = [('SID', '<i8', ()), ('V1', '<f8', ()), ('V2', '<f8', ()), ('ND', '<i8', ()), ('MSGLVL', '<i8', ()),
             ('MAXSET', '<i8', ()), ('SHFSCL', '<f8', ()), ('FLAG1', '<i8', ()), ('FLAG2', '<i8', ()),
             ('NORM', 'S4', ()), ('ALPH', '<f8', ()), ('FREQS_POS', '<i8', ()), ('FREQS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/EIGRL/FREQS']


@register_table
class FBADLAY(object):
    name = 'FBADLAY'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('PNTID', '<i8', ()), ('C', '<i8', ()),
             ('DELAY', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FBALOAD(object):
    name = 'FBALOAD'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('PNTID', '<i8', ()), ('C', '<i8', ()),
             ('A', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FBAPHAS(object):
    name = 'FBAPHAS'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('PNTID', '<i8', ()), ('C', '<i8', ()),
             ('PHASE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FREQS(object):
    name = 'FREQS'
    path = '/NASTRAN/INPUT/DYNAMIC/FREQ'
    dtype = [('F', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/FREQ'
    dtype = [('SID', '<i8', ()), ('FREQS_POS', '<i8', ()), ('FREQS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/FREQ/FREQS']


@register_table
class FREQ1(object):
    name = 'FREQ1'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('F1', '<f8', ()), ('DF', '<f8', ()), ('NDF', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FREQ2(object):
    name = 'FREQ2'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('NF', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FREQ3(object):
    name = 'FREQ3'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('TYPE', 'S4', ()), ('NEF', '<i8', ()),
             ('BIAS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FREQ4(object):
    name = 'FREQ4'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('TYPE', 'S4', ()), ('FSPD', '<f8', ()),
             ('NFM', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRACTIONS(object):
    name = 'FRACTIONS'
    path = '/NASTRAN/INPUT/DYNAMIC/FREQ5'
    dtype = [('FRI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/FREQ5'
    dtype = [('SID', '<i8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('FRACTIONS_POS', '<i8', ()),
             ('FRACTIONS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/FREQ5/FRACTIONS']


@register_table
class FRFCOMP(object):
    name = 'FRFCOMP'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('MEDIUM', '<i8', ()), ('UNITNO', '<i8', ()),
             ('LSCALEF', '<f8', ()), ('FSCALEF', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRFCONN(object):
    name = 'FRFCONN'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('COMPID1', '<i8', ()), ('COMPNAM1', 'S8', ()), ('PNTID1', '<i8', ()),
             ('COMPID2', '<i8', ()), ('COMPNAM2', 'S8', ()), ('PNTID2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRFFLEX(object):
    name = 'FRFFLEX'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('C', '<i8', ()), ('COMPID1', '<i8', ()), ('COMPNAM1', 'S8', ()),
             ('PNTID1', '<i8', ()), ('COMPID2', '<i8', ()), ('COMPNAM2', 'S8', ()), ('PNTID2', '<i8', ()),
             ('KVALUE', '<f8', ()), ('KTABID', '<i8', ()), ('BVALUE', '<f8', ()), ('BTABID', '<i8', ()),
             ('GEVALUE', '<f8', ()), ('GETABID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRFRELS(object):
    name = 'FRFRELS'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('C', '<i8', ()), ('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('GRIDID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRFSPC1(object):
    name = 'FRFSPC1'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('C', '<i8', ()), ('COMPID', '<i8', ()), ('COMPNAME', 'S8', ()), ('PNTID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRFXIT(object):
    name = 'FRFXIT'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('PNTID', '<i8', ()), ('C', '<i8', ()), ('LABEL', 'S48', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PNTLIST(object):
    name = 'PNTLIST'
    path = '/NASTRAN/INPUT/DYNAMIC/FRFXIT1'
    dtype = [('PNTID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/DYNAMIC/FRFXIT1'
    dtype = [('PNTID1', '<i8', ()), ('PNTID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/FRFXIT1'
    dtype = [('C', '<i8', ()), ('LIST_POS', '<i8', ()), ('LIST_LEN', '<i8', ()), ('THRU_POS', '<i8', ()),
             ('THRU_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/FRFXIT1/PNTLIST', '/NASTRAN/INPUT/DYNAMIC/FRFXIT1/THRU']


@register_table
class MONCARL(object):
    name = 'MONCARL'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('VARTYPE', 'S8', ()), ('STRVAR', 'S4', ()), ('EXISTSVR', 'S4', ()),
             ('MSVAR', '<f8', ()), ('KSVAR', '<f8', ()), ('BSVAR', '<f8', ()), ('FLUVAR', 'S4', ()),
             ('EXISTFVR', 'S4', ()), ('MFVAR', '<f8', ()), ('KFVAR', '<f8', ()), ('BFVAR', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NLRSFD(object):
    name = 'NLRSFD'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('PLANE', 'S8', ()), ('BDIA', '<f8', ()),
             ('BLEN', '<f8', ()), ('BCLR', '<f8', ()), ('SOLN', 'S8', ()), ('VISCO', '<f8', ()), ('PVAPCO', '<f8', ()),
             ('NPORT', '<i8', ()), ('PRES1', '<f8', ()), ('THETA1', '<f8', ()), ('PRES2', '<f8', ()),
             ('THETA2', '<f8', ()), ('NPNT', '<i8', ()), ('OFFSET1', '<f8', ()), ('OFFSET2', '<f8', ()),
             ('GRPNAME', 'S8', ()), ('NAME2', 'S8', ()), ('RDATA', '<f8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RANDPS(object):
    name = 'RANDPS'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('J', '<i8', ()), ('K', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()),
             ('TID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RANDT1(object):
    name = 'RANDT1'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('N', '<i8', ()), ('TO', '<f8', ()), ('TMAX', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RCROSS(object):
    name = 'RCROSS'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('RTYPE1', '<i8', ()), ('ID1', '<i8', ()), ('COMP1', '<i8', ()), ('RTYPE2', '<i8', ()),
             ('ID2', '<i8', ()), ('COMP2', '<i8', ()), ('BURID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RGYRO(object):
    name = 'RGYRO'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('ATYPE', 'S8', ()), ('REFROT', '<i8', ()), ('UNIT', 'S8', ()), ('SPDLOW', '<f8', ()),
             ('SPDHGH', '<f8', ()), ('SPEED', '<f8', ()), ('SEID', '<i8', ()), ('WR3WRL', '<f8', ()),
             ('WR4WRL', '<f8', ()), ('WRHWRL', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROTHYBD(object):
    name = 'ROTHYBD'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('ROTORID', '<i8', ()), ('HYBDMPID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ATHRU(object):
    name = 'ATHRU'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('GID1', '<i8', ()), ('GID2', '<i8', ()), ('BY', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class AXIS(object):
    name = 'AXIS'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('PID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class EIDS(object):
    name = 'EIDS'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('EID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class ETHRU(object):
    name = 'ETHRU'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('EID1', '<i8', ()), ('EID2', '<i8', ()), ('BY', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PIDS(object):
    name = 'PIDS'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('PID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PTHRU(object):
    name = 'PTHRU'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('PID1', '<i8', ()), ('PID2', '<i8', ()), ('BY', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTOR'
    dtype = [('ROTORID', '<i8', ()), ('FTYPE', 'S4', ()), ('EIDS_POS', '<i8', ()), ('EIDS_LEN', '<i8', ()),
             ('ETHRU_POS', '<i8', ()), ('ETHRU_LEN', '<i8', ()), ('PIDS_POS', '<i8', ()), ('PIDS_LEN', '<i8', ()),
             ('PTHRU_POS', '<i8', ()), ('PTHRU_LEN', '<i8', ()), ('AXIS_POS', '<i8', ()), ('AXIS_LEN', '<i8', ()),
             ('ATHRU_POS', '<i8', ()), ('ATHRU_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/ROTOR/ATHRU', '/NASTRAN/INPUT/DYNAMIC/ROTOR/AXIS',
                 '/NASTRAN/INPUT/DYNAMIC/ROTOR/EIDS', '/NASTRAN/INPUT/DYNAMIC/ROTOR/ETHRU',
                 '/NASTRAN/INPUT/DYNAMIC/ROTOR/PIDS', '/NASTRAN/INPUT/DYNAMIC/ROTOR/PTHRU']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTORAX'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTORAX'
    dtype = [('ROTORID', '<i8', ()), ('LISTTYPE', 'S4', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('INC', '<i8', ()),
             ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/ROTORAX/IDS']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTORG'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/ROTORG'
    dtype = [('ROTORID', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('INC', '<i8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/ROTORG/IDS']


@register_table
class ROTORSE(object):
    name = 'ROTORSE'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('RID', '<i8', ()), ('SEID', '<i8', ()), ('SEOPT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RSPINR(object):
    name = 'RSPINR'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('RID', '<i8', ()), ('GRIDA', '<i8', ()), ('GRIDB', '<i8', ()), ('UNIT', 'S8', ()), ('SPTID', '<i8', ()),
             ('SPTVAL', '<f8', ()), ('SEID', '<i8', ()), ('GR', '<f8', ()), ('ALPHAR1', '<f8', ()),
             ('ALPHAR2', '<f8', ()), ('WR3R', '<f8', ()), ('WR4R', '<f8', ()), ('WRHR', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RSPINT(object):
    name = 'RSPINT'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('RID', '<i8', ()), ('GRIDA', '<i8', ()), ('GRIDB', '<i8', ()), ('UNIT', 'S8', ()), ('TABLEID', '<i8', ()),
             ('SPDOUT', '<i8', ()), ('SEID', '<i8', ()), ('GR', '<f8', ()), ('ALPHAR1', '<f8', ()),
             ('ALPHAR2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TABLES(object):
    name = 'TABLES'
    path = '/NASTRAN/INPUT/DYNAMIC/SPECSEL'
    dtype = [('TID', '<i8', ()), ('DAMP', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/SPECSEL'
    dtype = [('ZERO', '<i8', ()), ('TYPE', 'S8', ()), ('TABLES_POS', '<i8', ()), ('TABLES_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/SPECSEL/TABLES']


@register_table
class SIGMA(object):
    name = 'SIGMA'
    path = '/NASTRAN/INPUT/DYNAMIC/TF'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('A0I', '<f8', ()), ('A1I', '<f8', ()), ('A2I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/TF'
    dtype = [('SID', '<i8', ()), ('GD', '<i8', ()), ('CD', '<i8', ()), ('B0', '<f8', ()), ('B1', '<f8', ()),
             ('B2', '<f8', ()), ('SIGMA_POS', '<i8', ()), ('SIGMA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/TF/SIGMA']


@register_table
class TIC(object):
    name = 'TIC'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('U0', '<f8', ()), ('V0', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/DYNAMIC/TMPSET'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/DYNAMIC/TMPSET'
    dtype = [('GROUPID', '<i8', ()), ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/DYNAMIC/TMPSET/GRIDS']


@register_table
class TTEMP(object):
    name = 'TTEMP'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('GROUPID', '<i8', ()), ('TID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class UNBALNC(object):
    name = 'UNBALNC'
    path = '/NASTRAN/INPUT/DYNAMIC'
    dtype = [('SID', '<i8', ()), ('MTABLE', '<i8', ()), ('MASS', '<f8', ()), ('GRID', '<i8', ()), ('X1', '<f8', ()),
             ('X2', '<f8', ()), ('X3', '<f8', ()), ('ROFFTAB', '<i8', ()), ('ROFFSET', '<f8', ()), ('THETA', '<f8', ()),
             ('ZOFFTAB', '<i8', ()), ('ZOFFSET', '<f8', ()), ('TON', '<f8', ()), ('TOFF', '<f8', ()),
             ('CFLAG', 'S8', ()), ('UFT1', '<i8', ()), ('UFT2', '<i8', ()), ('UFT3', '<i8', ()), ('UFR1', '<i8', ()),
             ('UFR2', '<i8', ()), ('UFR3', '<i8', ()), ('MCT1', '<i8', ()), ('MCT2', '<i8', ()), ('MCT3', '<i8', ()),
             ('MCR1', '<i8', ()), ('MCR2', '<i8', ()), ('MCR3', '<i8', ()), ('SCR1', '<i8', ()), ('SCR2', '<i8', ()),
             ('SCR3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEQUAD4(object):
    name = 'AEQUAD4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('CMPID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('G4', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEROD(object):
    name = 'AEROD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('AEID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('FIX', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEROQ4(object):
    name = 'AEROQ4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPTYPE', 'S8', ()), ('G', '<i8', (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEROT3(object):
    name = 'AEROT3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPTYPE', 'S8', ()), ('G', '<i8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AETRIA3(object):
    name = 'AETRIA3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('CMPID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAMAERO(object):
    name = 'BEAMAERO'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('COMPID', '<i8', ()), ('COMPTYPE', 'S8', ()), ('G', '<i8', (2,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BOTTOM(object):
    name = 'BOTTOM'
    path = '/NASTRAN/INPUT/ELEMENT/BOLT'
    dtype = [('GBID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TOP(object):
    name = 'TOP'
    path = '/NASTRAN/INPUT/ELEMENT/BOLT'
    dtype = [('GTID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/BOLT'
    dtype = [('ID', '<i8', ()), ('GRIDC', '<i8', ()), ('TOP_POS', '<i8', ()), ('TOP_LEN', '<i8', ()),
             ('BOTTOM_POS', '<i8', ()), ('BOTTOM_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/BOLT/BOTTOM', '/NASTRAN/INPUT/ELEMENT/BOLT/TOP']


@register_table
class CAABSF(object):
    name = 'CAABSF'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('G4', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CACINF3(object):
    name = 'CACINF3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CACINF4(object):
    name = 'CACINF4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('G4', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERO1(object):
    name = 'CAERO1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CP', '<i8', ()), ('NSPAN', '<i8', ()), ('NCHORD', '<i8', ()),
             ('LSPAN', '<i8', ()), ('LCHORD', '<i8', ()), ('IGID', '<i8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()),
             ('Z1', '<f8', ()), ('X12', '<f8', ()), ('X4', '<f8', ()), ('Y4', '<f8', ()), ('Z4', '<f8', ()),
             ('X43', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERO2(object):
    name = 'CAERO2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CP', '<i8', ()), ('NSB', '<i8', ()), ('NINT', '<i8', ()),
             ('LSB', '<i8', ()), ('LINT', '<i8', ()), ('IGID', '<i8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()),
             ('Z1', '<f8', ()), ('X12', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERO3(object):
    name = 'CAERO3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CP', '<i8', ()), ('LISTW', '<i8', ()), ('LISTC1', '<i8', ()),
             ('LISTC2', '<i8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('Z1', '<f8', ()), ('X12', '<f8', ()),
             ('X4', '<f8', ()), ('Y4', '<f8', ()), ('Z4', '<f8', ()), ('X43', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERO4(object):
    name = 'CAERO4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CP', '<i8', ()), ('NSPAN', '<i8', ()), ('LSPAN', '<i8', ()),
             ('X1', '<f8', ()), ('Y1', '<f8', ()), ('Z1', '<f8', ()), ('X12', '<f8', ()), ('X4', '<f8', ()),
             ('Y4', '<f8', ()), ('Z4', '<f8', ()), ('X43', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERO5(object):
    name = 'CAERO5'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CP', '<i8', ()), ('NSPAN', '<i8', ()), ('LSPAN', '<i8', ()),
             ('NTHRY', '<i8', ()), ('NTHICK', '<i8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('Z1', '<f8', ()),
             ('X12', '<f8', ()), ('X4', '<f8', ()), ('Y4', '<f8', ()), ('Z4', '<f8', ()), ('X43', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAXIF2(object):
    name = 'CAXIF2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('RHO', '<f8', ()), ('B', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAXIF3(object):
    name = 'CAXIF3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('IDF3', '<i8', ()), ('RHO', '<f8', ()),
             ('B', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAXIF4(object):
    name = 'CAXIF4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('IDF3', '<i8', ()), ('IDF4', '<i8', ()),
             ('RHO', '<f8', ()), ('B', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAXISYM(object):
    name = 'CAXISYM'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('OFF', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBAR(object):
    name = 'CBAR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('FLAG', '<i8', ()),
             ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('GO', '<i8', ()), ('PA', '<i8', ()),
             ('PB', '<i8', ()), ('W1A', '<f8', ()), ('W2A', '<f8', ()), ('W3A', '<f8', ()), ('W1B', '<f8', ()),
             ('W2B', '<f8', ()), ('W3B', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBEAM(object):
    name = 'CBEAM'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('SA', '<i8', ()),
             ('SB', '<i8', ()), ('X', '<f8', (3,)), ('G0', '<i8', ()), ('F', '<i8', ()), ('PA', '<i8', ()),
             ('PB', '<i8', ()), ('WA', '<f8', (3,)), ('WB', '<f8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBEAM3(object):
    name = 'CBEAM3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (3,)), ('DWRP', '<i8', (3,)), ('FLAG', '<i8', ()),
             ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('GO', '<i8', ()), ('W1A', '<f8', ()),
             ('W2A', '<f8', ()), ('W3A', '<f8', ()), ('W1B', '<f8', ()), ('W2B', '<f8', ()), ('W3B', '<f8', ()),
             ('W1C', '<f8', ()), ('W2C', '<f8', ()), ('W3C', '<f8', ()), ('TWA', '<f8', ()), ('TWB', '<f8', ()),
             ('TWC', '<f8', ()), ('PA', '<i8', ()), ('PB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBEND(object):
    name = 'CBEND'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('FLAG', '<i8', ()),
             ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('GO', '<i8', ()), ('GEOM', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBUSH(object):
    name = 'CBUSH'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('FLAG', '<i8', ()),
             ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('GO', '<i8', ()), ('CID', '<i8', ()),
             ('S', '<f8', ()), ('OCID', '<i8', ()), ('S1', '<f8', ()), ('S2', '<f8', ()), ('S3', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBUSH1D(object):
    name = 'CBUSH1D'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (2,)), ('CID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBUSH2D(object):
    name = 'CBUSH2D'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (2,)), ('CID', '<i8', ()), ('PLANE', '<i8', ()),
             ('SPTID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CCONEAX(object):
    name = 'CCONEAX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('RA', '<i8', ()), ('RB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDAMP1(object):
    name = 'CDAMP1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDAMP2(object):
    name = 'CDAMP2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('B', '<f8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDAMP3(object):
    name = 'CDAMP3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDAMP4(object):
    name = 'CDAMP4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('B', '<f8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CDAMP5(object):
    name = 'CDAMP5'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CELAS1(object):
    name = 'CELAS1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CELAS2(object):
    name = 'CELAS2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('K', '<f8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('GE', '<f8', ()), ('S', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CELAS3(object):
    name = 'CELAS3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CELAS4(object):
    name = 'CELAS4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('K', '<f8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFAST(object):
    name = 'CFAST'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GS', '<i8', ()), ('FORMAT', '<i8', ()), ('GA', '<i8', ()),
             ('GB', '<i8', ()), ('TYPE', '<i8', ()), ('CID', '<i8', ()), ('GUPPER', '<i8', (8,)),
             ('GLOWER', '<i8', (8,)), ('XYZ', '<f8', (3,)), ('TAVG', '<f8', ()), ('TMIN', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFASTP(object):
    name = 'CFASTP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('CID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()),
             ('TAVG', '<f8', ()), ('T1', '<f8', (3,)), ('T2', '<f8', (3,)), ('GH', '<i8', (8,)), ('GSHA', '<i8', ()),
             ('GSHB', '<i8', ()), ('GSHH', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFASTP_O(object):
    name = 'CFASTP_O'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('PRJTOL', '<f8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()),
             ('NMAP', '<i8', ()), ('CID', '<i8', ()), ('IECT', '<i8', (32,)), ('TAVG', '<f8', ()), ('TMIN', '<f8', ()),
             ('XYZS', '<f8', (3,)), ('XYZA', '<f8', (3,)), ('XYZB', '<f8', (3,)), ('NBIT', '<i8', (24,)),
             ('GIDA', '<i8', (4,)), ('GIDB', '<i8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFLUID2(object):
    name = 'CFLUID2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('RHO', '<f8', ()), ('B', '<f8', ()),
             ('HARMINDX', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFLUID3(object):
    name = 'CFLUID3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('IDF3', '<i8', ()), ('RHO', '<f8', ()),
             ('B', '<f8', ()), ('HARMINDX', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFLUID4(object):
    name = 'CFLUID4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDF1', '<i8', ()), ('IDF2', '<i8', ()), ('IDF3', '<i8', ()), ('IDF4', '<i8', ()),
             ('RHO', '<f8', ()), ('B', '<f8', ()), ('HARMINDX', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CGAP(object):
    name = 'CGAP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('F', '<i8', ()),
             ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('GO', '<i8', ()), ('CID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHACAB(object):
    name = 'CHACAB'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (20,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHACBR(object):
    name = 'CHACBR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (20,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHBDYE(object):
    name = 'CHBDYE'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('EID2', '<i8', ()), ('SIDE', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()),
             ('RADMIDF', '<i8', ()), ('RADMIDB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHBDYG(object):
    name = 'CHBDYG'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('TYPE', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()),
             ('RADMIDF', '<i8', ()), ('RADMIDB', '<i8', ()), ('G', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHBDYP(object):
    name = 'CHBDYP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('TYPE', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()),
             ('G1', '<i8', ()), ('G2', '<i8', ()), ('GO', '<i8', ()), ('RADMIDF', '<i8', ()), ('RADMIDB', '<i8', ()),
             ('DISLIN', '<i8', ()), ('CE', '<i8', ()), ('E1', '<f8', ()), ('E2', '<f8', ()), ('E3', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXA(object):
    name = 'CHEXA'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (20,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXAL(object):
    name = 'CHEXAL'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('MID', '<i8', ()), ('G', '<i8', (20,)), ('THETA', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXP(object):
    name = 'CHEXP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (8,)), ('E1', '<i8', (24,)), ('F', '<i8', (6,)),
             ('B1', '<i8', ()), ('E2', '<i8', (24,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CIFHEX(object):
    name = 'CIFHEX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (20,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CIFPENT(object):
    name = 'CIFPENT'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (15,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CIFQDX(object):
    name = 'CIFQDX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CIFQUAD(object):
    name = 'CIFQUAD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CMASS1(object):
    name = 'CMASS1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CMASS2(object):
    name = 'CMASS2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('M', '<f8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('C1', '<i8', ()),
             ('C2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CMASS3(object):
    name = 'CMASS3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CMASS4(object):
    name = 'CMASS4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('M', '<f8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONM2(object):
    name = 'CONM2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('G', '<i8', ()), ('CID', '<i8', ()), ('M', '<f8', ()), ('X1', '<f8', ()),
             ('X2', '<f8', ()), ('X3', '<f8', ()), ('I1', '<f8', ()), ('I2', '<f8', (2,)), ('I3', '<f8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD(object):
    name = 'CONROD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('MID', '<i8', ()), ('A', '<f8', ()),
             ('J', '<f8', ()), ('C', '<f8', ()), ('NSM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONTRLT(object):
    name = 'CONTRLT'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('CNTLID', '<i8', ()), ('SENSOR', '<i8', ()), ('SFORM', '<i8', ()), ('CTYPE', '<i8', ()),
             ('PLOW', '<f8', ()), ('PHIGH', '<f8', ()), ('PTYPE', '<i8', ()), ('PZERO', '<f8', ()), ('DT', '<f8', ()),
             ('DELAY', '<f8', ()), ('TAUD', '<f8', ()), ('TAUI', '<f8', ()), ('SETPT', '<f8', ()), ('ITID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CPENTA(object):
    name = 'CPENTA'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (15,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQDX4FD(object):
    name = 'CQDX4FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (9,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQDX9FD(object):
    name = 'CQDX9FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (9,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD(object):
    name = 'CQUAD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (9,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4(object):
    name = 'CQUAD4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (4,)), ('THETA', '<f8', ()), ('ZOFFS', '<f8', ()),
             ('TFLAG', '<i8', ()), ('T', '<f8', (4,)), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4FD(object):
    name = 'CQUAD4FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (9,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD8(object):
    name = 'CQUAD8'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (8,)), ('T', '<f8', (4,)), ('THETA', '<f8', ()),
             ('ZOFFS', '<f8', ()), ('TFLAG', '<i8', ()), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD9FD(object):
    name = 'CQUAD9FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (9,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUADR(object):
    name = 'CQUADR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (4,)), ('THETA', '<f8', ()), ('ZOFFS', '<f8', ()),
             ('TFLAG', '<i8', ()), ('T', '<f8', (4,)), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUADX(object):
    name = 'CQUADX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/INPUT/ELEMENT/CQUAD'
    subtables = []


@register_table
class DEPEN(object):
    name = 'DEPEN'
    path = '/NASTRAN/INPUT/ELEMENT/CRBE1'
    dtype = [('GM', '<i8', ()), ('CM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class INDEP(object):
    name = 'INDEP'
    path = '/NASTRAN/INPUT/ELEMENT/CRBE1'
    dtype = [('GN', '<i8', ()), ('CN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class LAGMULTIPL(object):
    name = 'LAGMULTIPL'
    path = '/NASTRAN/INPUT/ELEMENT/CRBE1'
    dtype = [('LMID', '<i8', ()), ('NDOF', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THERMEXP(object):
    name = 'THERMEXP'
    path = '/NASTRAN/INPUT/ELEMENT/CRBE1'
    dtype = [('ALPHA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/CRBE1'
    dtype = [('EID', '<i8', ()), ('NWE', '<i8', ()), ('ELTYPE', '<i8', ()), ('INDEP_POS', '<i8', ()),
             ('INDEP_LEN', '<i8', ()), ('DEPEN_POS', '<i8', ()), ('DEPEN_LEN', '<i8', ()), ('THERMEXP_POS', '<i8', ()),
             ('THERMEXP_LEN', '<i8', ()), ('LAGM_POS', '<i8', ()), ('LAGM_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/CRBE1/DEPEN', '/NASTRAN/INPUT/ELEMENT/CRBE1/INDEP',
                 '/NASTRAN/INPUT/ELEMENT/CRBE1/LAGMULTIPL', '/NASTRAN/INPUT/ELEMENT/CRBE1/THERMEXP']


@register_table
class CROD(object):
    name = 'CROD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CSEAM(object):
    name = 'CSEAM'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('SMLN', 'S8', ()), ('SMLNID', '<i8', (2,)), ('CTYPE', '<i8', ()),
             ('IDAS', '<i8', ()), ('IDBS', '<i8', ()), ('IDAE', '<i8', ()), ('IDBE', '<i8', ()), ('GFLAG', '<i8', ()),
             ('GS', '<i8', ()), ('GE', '<i8', ()), ('XS', '<f8', ()), ('YS', '<f8', ()), ('ZS', '<f8', ()),
             ('XE', '<f8', ()), ('YE', '<f8', ()), ('ZE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CSEAMP(object):
    name = 'CSEAMP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('SMLN', 'S8', ()), ('GS', '<i8', ()), ('GE', '<i8', ()),
             ('GAX', '<i8', (8,)), ('XYZS', '<f8', (3,)), ('XYZE', '<f8', (3,)), ('ISHR', '<i8', ()),
             ('GSHH', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CSHEAR(object):
    name = 'CSHEAR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CSLOT3(object):
    name = 'CSLOT3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDS', '<i8', (3,)), ('RHO', '<f8', ()), ('B', '<f8', ()), ('M', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CSLOT4(object):
    name = 'CSLOT4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('IDS', '<i8', (4,)), ('RHO', '<f8', ()), ('B', '<f8', ()), ('M', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTETRA(object):
    name = 'CTETRA'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (10,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA3(object):
    name = 'CTRIA3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (3,)), ('THETA', '<f8', ()), ('ZOFFS', '<f8', ()),
             ('TFLAG', '<i8', ()), ('T', '<f8', (3,)), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA3FD(object):
    name = 'CTRIA3FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA6(object):
    name = 'CTRIA6'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('ZOFFS', '<f8', ()),
             ('T', '<f8', (3,)), ('TFLAG', '<i8', ()), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA6FD(object):
    name = 'CTRIA6FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAH(object):
    name = 'CTRIAH'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/INPUT/ELEMENT/CTRIAX'
    subtables = []


@register_table
class CTRIAR(object):
    name = 'CTRIAR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (3,)), ('THETA', '<f8', ()), ('ZOFFS', '<f8', ()),
             ('TFLAG', '<i8', ()), ('T', '<f8', (3,)), ('MCID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAX(object):
    name = 'CTRIAX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAX6(object):
    name = 'CTRIAX6'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('MID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIX3FD(object):
    name = 'CTRIX3FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIX6FD(object):
    name = 'CTRIX6FD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (6,)), ('THETA', '<f8', ()), ('MCID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTUBE(object):
    name = 'CTUBE'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CVISC(object):
    name = 'CVISC'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('G', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CWELD(object):
    name = 'CWELD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GS', '<i8', ()), ('FORMAT', '<i8', ()), ('GA', '<i8', ()),
             ('GB', '<i8', ()), ('TYPE', '<i8', ()), ('MCID', '<i8', ()), ('GUPPER', '<i8', (8,)),
             ('GLOWER', '<i8', (8,)), ('SHIDA', '<i8', ()), ('SHIDB', '<i8', ()), ('XYZ', '<f8', (3,)),
             ('PIDA', '<i8', ()), ('PIDB', '<i8', ()), ('TAVG', '<f8', ()), ('TMIN', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CWELDC(object):
    name = 'CWELDC'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GS', '<i8', ()), ('FORMAT', '<i8', ()), ('GA', '<i8', ()),
             ('GB', '<i8', ()), ('TYPE', '<i8', ()), ('MCID', '<i8', ()), ('GUPPER', '<i8', (8,)),
             ('GLOWER', '<i8', (8,)), ('TAVG', '<f8', ()), ('PRJTOL', '<f8', ()), ('TMIN', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CWELDP(object):
    name = 'CWELDP'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('MCID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()),
             ('TAVG', '<f8', ()), ('TMIN', '<f8', ()), ('T3', '<f8', (3,)), ('GH', '<i8', (8,)), ('GSHA', '<i8', ()),
             ('GSHB', '<i8', ()), ('GSHH', '<i8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CWELDP_O(object):
    name = 'CWELDP_O'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('PRJTOL', '<f8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()),
             ('NMAP', '<i8', ()), ('MCID', '<i8', ()), ('IECT', '<i8', (32,)), ('TAVG', '<f8', ()), ('TMIN', '<f8', ()),
             ('XYZS', '<f8', (3,)), ('XYZA', '<f8', (3,)), ('XYZB', '<f8', (3,)), ('NBIT', '<i8', (24,)),
             ('GIDA', '<i8', (4,)), ('GIDB', '<i8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class KORZ(object):
    name = 'KORZ'
    path = '/NASTRAN/INPUT/ELEMENT/GENEL'
    dtype = [('KZIJ', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class S(object):
    name = 'S'
    path = '/NASTRAN/INPUT/ELEMENT/GENEL'
    dtype = [('SIJ', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class UDLIST(object):
    name = 'UDLIST'
    path = '/NASTRAN/INPUT/ELEMENT/GENEL'
    dtype = [('UD', '<i8', ()), ('CD', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class UILIST(object):
    name = 'UILIST'
    path = '/NASTRAN/INPUT/ELEMENT/GENEL'
    dtype = [('UI', '<i8', ()), ('CI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/GENEL'
    dtype = [('EID', '<i8', ()), ('M', '<i8', ()), ('N', '<i8', ()), ('F', '<i8', ()), ('ZEROS', '<i8', ()),
             ('UI_POS', '<i8', ()), ('UI_LEN', '<i8', ()), ('UD_POS', '<i8', ()), ('UD_LEN', '<i8', ()),
             ('NLT', '<i8', ()), ('KORZ_POS', '<i8', ()), ('S_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/GENEL/KORZ', '/NASTRAN/INPUT/ELEMENT/GENEL/S',
                 '/NASTRAN/INPUT/ELEMENT/GENEL/UDLIST', '/NASTRAN/INPUT/ELEMENT/GENEL/UILIST']


@register_table
class PLOTEL(object):
    name = 'PLOTEL'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('G', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM1(object):
    name = 'PRIM1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM2(object):
    name = 'PRIM2'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('P41', '<f8', ()), ('P42', '<f8', ()), ('P43', '<f8', ()), ('AMESH', '<i8', ()),
             ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM3(object):
    name = 'PRIM3'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM4(object):
    name = 'PRIM4'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('DIAM1', '<f8', ()), ('DIAM2', '<f8', ()), ('ANGLE1', '<f8', ()),
             ('ANGLE2', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM5(object):
    name = 'PRIM5'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('DIAM1', '<f8', ()), ('ANGLE1', '<f8', ()), ('ANGLE2', '<f8', ()),
             ('AMESH', '<i8', ()), ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM6(object):
    name = 'PRIM6'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('DIAM1', '<f8', ()), ('DIAM2', '<f8', ()), ('ANGLE1', '<f8', ()),
             ('ANGLE2', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM7(object):
    name = 'PRIM7'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('DIAM1', '<f8', ()), ('ANGLE1', '<f8', ()), ('ANGLE2', '<f8', ()),
             ('TRUNC1', '<f8', ()), ('TRUNC2', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIM8(object):
    name = 'PRIM8'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('PRIMID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('P11', '<f8', ()), ('P12', '<f8', ()), ('P13', '<f8', ()),
             ('P21', '<f8', ()), ('P22', '<f8', ()), ('P23', '<f8', ()), ('P31', '<f8', ()), ('P32', '<f8', ()),
             ('P33', '<f8', ()), ('DIAM1', '<f8', ()), ('ANGLE1', '<f8', ()), ('ANGLE2', '<f8', ()),
             ('TRUNC1', '<f8', ()), ('TRUNC2', '<f8', ()), ('AMESH', '<i8', ()), ('BMESH', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RADCOL(object):
    name = 'RADCOL'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('RCOLID', '<i8', ()), ('IVIEWF', '<i8', ()), ('IVIEWB', '<i8', ()), ('RADMIDF', '<i8', ()),
             ('RADMIDB', '<i8', ()), ('SET3ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RBAR(object):
    name = 'RBAR'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('CNA', '<i8', ()), ('CNB', '<i8', ()),
             ('CMA', '<i8', ()), ('CMB', '<i8', ()), ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RBAR1(object):
    name = 'RBAR1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('CMB', '<i8', ()), ('ALPHA', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DEP(object):
    name = 'DEP'
    path = '/NASTRAN/INPUT/ELEMENT/RBE1'
    dtype = [('GM', '<i8', ()), ('CM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class INDEP(object):
    name = 'INDEP'
    path = '/NASTRAN/INPUT/ELEMENT/RBE1'
    dtype = [('GN', '<i8', ()), ('CN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/RBE1'
    dtype = [('EID', '<i8', ()), ('ALPHA', '<f8', ()), ('INDEP_LEN', '<i8', ()), ('INDEP_POS', '<i8', ()),
             ('DEP_LEN', '<i8', ()), ('DEP_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/RBE1/DEP', '/NASTRAN/INPUT/ELEMENT/RBE1/INDEP']


@register_table
class GM(object):
    name = 'GM'
    path = '/NASTRAN/INPUT/ELEMENT/RBE2'
    dtype = [('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RB(object):
    name = 'RB'
    path = '/NASTRAN/INPUT/ELEMENT/RBE2'
    dtype = [('EID', '<i8', ()), ('GN', '<i8', ()), ('CM', '<i8', ()), ('GM_POS', '<i8', ()), ('GM_LEN', '<i8', ()),
             ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GMS(object):
    name = 'GMS'
    path = '/NASTRAN/INPUT/ELEMENT/RBE2GS'
    dtype = [('GM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GNS(object):
    name = 'GNS'
    path = '/NASTRAN/INPUT/ELEMENT/RBE2GS'
    dtype = [('GN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/RBE2GS'
    dtype = [('EID', '<i8', ()), ('GS', '<i8', ()), ('TYPE', '<i8', ()), ('IDP', '<i8', ()), ('SEID', '<i8', ()),
             ('CMB', '<i8', ()), ('R', '<f8', ()), ('ALPHA', '<f8', ()), ('XS', '<f8', ()), ('YS', '<f8', ()),
             ('ZS', '<f8', ()), ('GNS_POS', '<i8', ()), ('GNS_LEN', '<i8', ()), ('GMS_POS', '<i8', ()),
             ('GMS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/RBE2GS/GMS', '/NASTRAN/INPUT/ELEMENT/RBE2GS/GNS']


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/ELEMENT/RBE3'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GM(object):
    name = 'GM'
    path = '/NASTRAN/INPUT/ELEMENT/RBE3'
    dtype = [('GM', '<i8', ()), ('CM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class WTCG(object):
    name = 'WTCG'
    path = '/NASTRAN/INPUT/ELEMENT/RBE3'
    dtype = [('WT1', '<f8', ()), ('C', '<i8', ()), ('G_POS', '<i8', ()), ('G_LEN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/RBE3'
    dtype = [('EID', '<i8', ()), ('REFG', '<i8', ()), ('REFC', '<i8', ()), ('WTCG_POS', '<i8', ()),
             ('WTCG_LEN', '<i8', ()), ('GM_POS', '<i8', ()), ('GM_LEN', '<i8', ()), ('ALPHA', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/RBE3/G', '/NASTRAN/INPUT/ELEMENT/RBE3/GM', '/NASTRAN/INPUT/ELEMENT/RBE3/WTCG']


@register_table
class RINGAX(object):
    name = 'RINGAX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('ID', '<i8', ()), ('R', '<f8', ()), ('Z', '<f8', ()), ('PS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RJOINT(object):
    name = 'RJOINT'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('CB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RROD(object):
    name = 'RROD'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('CMA', '<i8', ()), ('CMB', '<i8', ()),
             ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RSSCON(object):
    name = 'RSSCON'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('STYLE', '<i8', ()), ('GRID1', '<i8', ()), ('GRID2', '<i8', ()), ('GRID3', '<i8', ()),
             ('GRID4', '<i8', ()), ('GRID5', '<i8', ()), ('GRID6', '<i8', ()), ('EDGE1', '<i8', ()),
             ('EDGE2', '<i8', ()), ('EDGE3', '<i8', ()), ('ELID1', '<i8', ()), ('ELID2', '<i8', ()),
             ('CBID', '<i8', ()), ('SBID', '<i8', ()), ('CBPID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTRPLT(object):
    name = 'RTRPLT'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('GC', '<i8', ()), ('CNA', '<i8', ()),
             ('CNB', '<i8', ()), ('CNC', '<i8', ()), ('CMA', '<i8', ()), ('CMB', '<i8', ()), ('CMC', '<i8', ()),
             ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTRPLT1(object):
    name = 'RTRPLT1'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('EID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('GC', '<i8', ()), ('CMB', '<i8', ()),
             ('CMC', '<i8', ()), ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTAX(object):
    name = 'SECTAX'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('ID', '<i8', ()), ('RID', '<i8', ()), ('R', '<f8', ()), ('PHI1', '<f8', ()), ('PHI2', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WETELME(object):
    name = 'WETELME'
    path = '/NASTRAN/INPUT/ELEMENT'
    dtype = [('WEID', '<i8', ()), ('EID', '<i8', ()), ('SIDE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GIDS(object):
    name = 'GIDS'
    path = '/NASTRAN/INPUT/ELEMENT/WETELMG'
    dtype = [('GID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/ELEMENT/WETELMG'
    dtype = [('WEID', '<i8', ()), ('TYPE', 'S8', ()), ('GIDS_LEN', '<i8', ()), ('GIDS_POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/ELEMENT/WETELMG/GIDS']


@register_table
class ELSET(object):
    name = 'ELSET'
    path = '/NASTRAN/INPUT/FATIGUE/FTGDEF'
    dtype = [('ELSID', '<i8', ()), ('PFTGID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SEAMW(object):
    name = 'SEAMW'
    path = '/NASTRAN/INPUT/FATIGUE/FTGDEF'
    dtype = [('ELSID', '<i8', ()), ('PFTGID', '<i8', ()), ('NDSID', '<i8', ()), ('WELD', 'S8', ()), ('TYPE', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SPOTW(object):
    name = 'SPOTW'
    path = '/NASTRAN/INPUT/FATIGUE/FTGDEF'
    dtype = [('ELSID', '<i8', ()), ('PFTGID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class XELSET(object):
    name = 'XELSET'
    path = '/NASTRAN/INPUT/FATIGUE/FTGDEF'
    dtype = [('XELSID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/FATIGUE/FTGDEF'
    dtype = [('ID', '<i8', ()), ('TOPSTR', '<f8', ()), ('PFTGID', '<i8', ()), ('TOPDMG', '<f8', ()),
             ('ELSET_POS', '<i8', ()), ('ELSET_LEN', '<i8', ()), ('SPOTW_POS', '<i8', ()), ('SPOTW_LEN', '<i8', ()),
             ('SEAMW_POS', '<i8', ()), ('SEAMW_LEN', '<i8', ()), ('XELSET_POS', '<i8', ()), ('XELSET_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/FATIGUE/FTGDEF/ELSET', '/NASTRAN/INPUT/FATIGUE/FTGDEF/SEAMW',
                 '/NASTRAN/INPUT/FATIGUE/FTGDEF/SPOTW', '/NASTRAN/INPUT/FATIGUE/FTGDEF/XELSET']


@register_table
class FLOAD(object):
    name = 'FLOAD'
    path = '/NASTRAN/INPUT/FATIGUE/FTGEVNT'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/FATIGUE/FTGEVNT'
    dtype = [('ID', '<i8', ()), ('FLOAD_POS', '<i8', ()), ('FLOAD_LEN', '<i8', ()), ('NAME', 'S4', ()),
             ('EXIST', 'S4', ()), ('EQNAME', 'S56', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/FATIGUE/FTGEVNT/FLOAD']


@register_table
class FTGLOAD(object):
    name = 'FTGLOAD'
    path = '/NASTRAN/INPUT/FATIGUE'
    dtype = [('ID', '<i8', ()), ('TID', '<i8', ()), ('LCID', '<i8', ()), ('LDM', '<f8', ()), ('SCALE_MAX', '<f8', ()),
             ('OFFSET_MIN', '<f8', ()), ('TYPE', 'S8', ()), ('CHNL', '<i8', ()), ('UNITS', 'S4', ()),
             ('EXIST', 'S4', ()), ('EQUIV', '<f8', ()), ('EQNAME', 'S48', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FTGPARM(object):
    name = 'FTGPARM'
    path = '/NASTRAN/INPUT/FATIGUE'
    dtype = [('ID', '<i8', ()), ('TYPE', 'S4', ()), ('FACTOR', '<f8', ()), ('NTHRD', '<i8', ()), ('LOGLVL', '<i8', ()),
             ('LAYER', '<i8', ()), ('STROUT', '<i8', ()), ('STRESS', 'S4', ()), ('existSTS', 'S4', ()),
             ('COMBS', 'S8', ()), ('CORRS', 'S8', ()), ('PLASTS', 'S8', ()), ('LOCS', 'S4', ()), ('INTERPS', '<i8', ()),
             ('RECS', 'S8', ()), ('SRESS', 'S4', ()), ('STRAIN', 'S4', ()), ('existSTN', 'S4', ()), ('COMBN', 'S8', ()),
             ('CORRN', 'S8', ()), ('PLASTN', 'S8', ()), ('LOCN', 'S4', ()), ('INTERPN', '<i8', ()), ('RECN', 'S8', ()),
             ('SRESN', 'S4', ()), ('RAINFLOW', 'S4', ()), ('existRNF', 'S4', ()), ('RTYPE', 'S8', ()),
             ('GATE', '<f8', ()), ('PCTRED', '<f8', ()), ('CERTNTY', 'S4', ()), ('existCRT', 'S4', ()),
             ('SURV', '<f8', ()), ('FOS', 'S4', ()), ('existFOS', 'S4', ()), ('OPTION', 'S8', ()), ('LIFE', '<f8', ()),
             ('BACKACC', '<f8', ()), ('MAXFAC', '<f8', ()), ('MINFAC', '<f8', ()), ('DAMAGE', 'S4', ()),
             ('existDAM', 'S4', ()), ('CHECK', '<i8', ()), ('FLOOR', '<f8', ()), ('MAXDAM', '<f8', ()),
             ('FAILDAM', '<f8', ()), ('SPOTW', 'S4', ()), ('existSPW', 'S4', ()), ('COMBW', 'S8', ()),
             ('CORRW', 'S8', ()), ('NANGLE', '<i8', ()), ('SWLOC', '<i8', ()), ('MIDDLE', '<i8', ()),
             ('TORSION', '<i8', ()), ('SEAMW', 'S4', ()), ('existSEW', 'S4', ()), ('COMBE', 'S8', ()),
             ('CORRE', 'S8', ()), ('THICK', '<i8', ()), ('LOCSM', 'S8', ()), ('RESENT', 'S8', ()), ('MULTI', 'S4', ()),
             ('existMLT', 'S4', ()), ('MMTHD', 'S8', ()), ('NONLWR', '<f8', ()), ('NONUPR', '<f8', ()),
             ('BIAXLWR', '<f8', ()), ('BIAXMID', '<f8', ()), ('BIAXUPR', '<f8', ()), ('ZERO', '<f8', ()),
             ('GATEM', '<f8', ()), ('NAVG', 'S4', ()), ('existNAV', 'S4', ()), ('MTHD', 'S4', ()),
             ('OUTPUT', '<i8', ()), ('NORMAL', 'S4', ()), ('VIBFTG', 'S4', ()), ('existVIB', 'S4', ()),
             ('ATYPE', 'S8', ()), ('MAXSTR', '<f8', ()), ('CLIPLVL', '<f8', ()), ('MAXPEAK', '<f8', ()),
             ('STRBINS', '<i8', ()), ('MAXFREQ', '<f8', ()), ('NCALC', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FID_N(object):
    name = 'FID_N'
    path = '/NASTRAN/INPUT/FATIGUE/FTGSEQ'
    dtype = [('FID', '<i8', ()), ('N', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/FATIGUE/FTGSEQ'
    dtype = [('ID', '<i8', ()), ('EVNTOUT', '<i8', ()), ('METHOD', '<i8', ()), ('TUNIT', 'S4', ()), ('LDM', '<f8', ()),
             ('SCALE', '<f8', ()), ('OFFSET', '<f8', ()), ('FID_N_POS', '<i8', ()), ('FID_N_LEN', '<i8', ()),
             ('UNITS', 'S4', ()), ('EXIST', 'S4', ()), ('EQUIV', '<f8', ()), ('EQNAME', 'S48', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/FATIGUE/FTGSEQ/FID_N']


@register_table
class MATID(object):
    name = 'MATID'
    path = '/NASTRAN/INPUT/FATIGUE/MATFTG'
    dtype = [('MID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TABLE(object):
    name = 'TABLE'
    path = '/NASTRAN/INPUT/FATIGUE/MATFTG'
    dtype = [('VALUE', '<f8', ()), ('TID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/FATIGUE/MATFTG'
    dtype = [('ID', '<i8', ()), ('CNVRT', '<f8', ()), ('STATIC', 'S4', ()), ('existSTA', 'S4', ()), ('YS', '<f8', ()),
             ('UTS', '<f8', ()), ('CODE', '<i8', ()), ('TYPE', 'S8', ()), ('RR', '<f8', ()), ('SE', '<f8', ()),
             ('MP', '<f8', ()), ('RA', '<f8', ()), ('SN', 'S4', ()), ('existSN', 'S4', ()), ('SRI1', '<f8', ()),
             ('b1', '<f8', ()), ('Nc1', '<f8', ()), ('b2', '<f8', ()), ('Nfc', '<f8', ()), ('SE_SN', '<f8', ()),
             ('BTHRESH', '<f8', ()), ('M1', '<f8', ()), ('M2', '<f8', ()), ('M3', '<f8', ()), ('M4', '<f8', ()),
             ('MSS', '<f8', ()), ('RTHICK', '<f8', ()), ('NTHICK', '<f8', ()), ('SF_FXY', '<f8', ()),
             ('DE_FXY', '<f8', ()), ('TE_FXY', '<f8', ()), ('SF_MXY', '<f8', ()), ('DE_MXY', '<f8', ()),
             ('TE_MXY', '<f8', ()), ('SF_FZ', '<f8', ()), ('DE_FZ', '<f8', ()), ('TE_FZ', '<f8', ()),
             ('SF_MZ', '<f8', ()), ('DE_MZ', '<f8', ()), ('TE_MZ', '<f8', ()), ('SNS1', 'S4', ()),
             ('exisSNS1', 'S4', ()), ('SRI1_S1', '<f8', ()), ('b1_S1', '<f8', ()), ('Nc1_S1', '<f8', ()),
             ('b2_S1', '<f8', ()), ('Nfc_S1', '<f8', ()), ('SE_SN_S1', '<f8', ()), ('BTHRESH_S1', '<f8', ()),
             ('M1_S1', '<f8', ()), ('M2_S1', '<f8', ()), ('M3_S1', '<f8', ()), ('M4_S1', '<f8', ()),
             ('MSS_S1', '<f8', ()), ('RTHICK_S1', '<f8', ()), ('NTHICK_S1', '<f8', ()), ('SF_FXY_S1', '<f8', ()),
             ('DE_FXY_S1', '<f8', ()), ('TE_FXY_S1', '<f8', ()), ('SF_MXY_S1', '<f8', ()), ('DE_MXY_S1', '<f8', ()),
             ('TE_MXY_S1', '<f8', ()), ('SF_FZ_S1', '<f8', ()), ('DE_FZ_S1', '<f8', ()), ('TE_FZ_S1', '<f8', ()),
             ('SF_MZ_S1', '<f8', ()), ('DE_MZ_S1', '<f8', ()), ('TE_MZ_S1', '<f8', ()), ('SNS2', 'S4', ()),
             ('exisSNS2', 'S4', ()), ('SRI1_S2', '<f8', ()), ('b1_S2', '<f8', ()), ('Nc1_S2', '<f8', ()),
             ('b2_S2', '<f8', ()), ('Nfc_S2', '<f8', ()), ('SE_SN_S2', '<f8', ()), ('BTHRESH_S2', '<f8', ()),
             ('M1_S2', '<f8', ()), ('M2_S2', '<f8', ()), ('M3_S2', '<f8', ()), ('M4_S2', '<f8', ()),
             ('MSS_S2', '<f8', ()), ('RTHICK_S2', '<f8', ()), ('NTHICK_S2', '<f8', ()), ('SF_FXY_S2', '<f8', ()),
             ('DE_FXY_S2', '<f8', ()), ('TE_FXY_S2', '<f8', ()), ('SF_MXY_S2', '<f8', ()), ('DE_MXY_S2', '<f8', ()),
             ('TE_MXY_S2', '<f8', ()), ('SF_FZ_S2', '<f8', ()), ('DE_FZ_S2', '<f8', ()), ('TE_FZ_S2', '<f8', ()),
             ('SF_MZ_S2', '<f8', ()), ('DE_MZ_S2', '<f8', ()), ('TE_MZ_S2', '<f8', ()), ('SNBR0', 'S4', ()),
             ('exisSNBR0', 'S4', ()), ('SRI1_R0', '<f8', ()), ('b1_R0', '<f8', ()), ('Nc1_R0', '<f8', ()),
             ('b2_R0', '<f8', ()), ('Nfc_R0', '<f8', ()), ('SE_SN_R0', '<f8', ()), ('BTHRESH_R0', '<f8', ()),
             ('M1_R0', '<f8', ()), ('M2_R0', '<f8', ()), ('M3_R0', '<f8', ()), ('M4_R0', '<f8', ()),
             ('MSS_R0', '<f8', ()), ('RTHICK_R0', '<f8', ()), ('NTHICK_R0', '<f8', ()), ('SF_FXY_R0', '<f8', ()),
             ('DE_FXY_R0', '<f8', ()), ('TE_FXY_R0', '<f8', ()), ('SF_MXY_R0', '<f8', ()), ('DE_MXY_R0', '<f8', ()),
             ('TE_MXY_R0', '<f8', ()), ('SF_FZ_R0', '<f8', ()), ('DE_FZ_R0', '<f8', ()), ('TE_FZ_R0', '<f8', ()),
             ('SF_MZ_R0', '<f8', ()), ('DE_MZ_R0', '<f8', ()), ('TE_MZ_R0', '<f8', ()), ('SNBR1', 'S4', ()),
             ('exisSNBR1', 'S4', ()), ('SRI1_R1', '<f8', ()), ('b1_R1', '<f8', ()), ('Nc1_R1', '<f8', ()),
             ('b2_R1', '<f8', ()), ('Nfc_R1', '<f8', ()), ('SE_SN_R1', '<f8', ()), ('BTHRESH_R1', '<f8', ()),
             ('M1_R1', '<f8', ()), ('M2_R1', '<f8', ()), ('M3_R1', '<f8', ()), ('M4_R1', '<f8', ()),
             ('MSS_R1', '<f8', ()), ('RTHICK_R1', '<f8', ()), ('NTHICK_R1', '<f8', ()), ('SF_FXY_R1', '<f8', ()),
             ('DE_FXY_R1', '<f8', ()), ('TE_FXY_R1', '<f8', ()), ('SF_MXY_R1', '<f8', ()), ('DE_MXY_R1', '<f8', ()),
             ('TE_MXY_R1', '<f8', ()), ('SF_FZ_R1', '<f8', ()), ('DE_FZ_R1', '<f8', ()), ('TE_FZ_R1', '<f8', ()),
             ('SF_MZ_R1', '<f8', ()), ('DE_MZ_R1', '<f8', ()), ('TE_MZ_R1', '<f8', ()), ('TABL', 'S4', ()),
             ('TABLE_POS', '<i8', ()), ('TABLE_LEN', '<i8', ()), ('BASTEN', 'S4', ()), ('existBAS', 'S4', ()),
             ('BT_A', '<f8', ()), ('BT_B', '<f8', ()), ('BT_c', '<f8', ()), ('BT_Eb', '<f8', ()), ('BT_Sd', '<f8', ()),
             ('EN', 'S4', ()), ('existEN', 'S4', ()), ('EN_Sf', '<f8', ()), ('EN_b', '<f8', ()), ('EN_c', '<f8', ()),
             ('EN_Ef', '<f8', ()), ('EN_n', '<f8', ()), ('EN_K', '<f8', ()), ('EN_Nc', '<f8', ()),
             ('EN_SEe', '<f8', ()), ('EN_SEp', '<f8', ()), ('EN_SEc', '<f8', ()), ('EN_Ne', '<f8', ()),
             ('EN_FSN', '<f8', ()), ('EN_S', '<f8', ()), ('MATID', 'S4', ()), ('E', '<f8', ()), ('NU', '<f8', ()),
             ('MATID_POS', '<i8', ()), ('MATID_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/FATIGUE/MATFTG/MATID', '/NASTRAN/INPUT/FATIGUE/MATFTG/TABLE']


@register_table
class PFTG(object):
    name = 'PFTG'
    path = '/NASTRAN/INPUT/FATIGUE'
    dtype = [('ID', '<i8', ()), ('LAYER', '<i8', ()), ('FINISH', 'S8', ()), ('KFINISH', '<f8', ()), ('KF', '<f8', ()),
             ('SCALE', '<f8', ()), ('OFFSET', '<f8', ()), ('SHAPE', '<f8', ()), ('KTREAT', '<f8', ()),
             ('DIAM', '<f8', ()), ('T1', '<f8', ()), ('T2', '<f8', ()), ('SPTFLG', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EVNT(object):
    name = 'EVNT'
    path = '/NASTRAN/INPUT/FATIGUE/TIM2PSD'
    dtype = [('EVNTID', '<i8', ()), ('NSI', '<i8', ()), ('RMSI', '<i8', ()), ('TSMOOTH', '<i8', ()), ('SF', '<f8', ()),
             ('T', '<f8', ()), ('DELTA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DEL(object):
    name = 'DEL'
    path = '/NASTRAN/INPUT/FATIGUE/TIM2PSD'
    dtype = [('EVENID', '<i8', ()), ('TI1', '<f8', ()), ('TF1', '<f8', ()), ('TI2', '<f8', ()), ('TF2', '<f8', ()),
             ('TI3', '<f8', ()), ('TF3', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class MAP(object):
    name = 'MAP'
    path = '/NASTRAN/INPUT/FATIGUE/TIM2PSD'
    dtype = [('LCIDI', '<i8', ()), ('CHANI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/FATIGUE/TIM2PSD'
    dtype = [('ID', '<i8', ()), ('SRATE', '<f8', ()), ('WINDOW', 'S8', ()), ('FORMAT', 'S4', ()), ('MEANS', 'S4', ()),
             ('NSKIP', '<i8', ()), ('MAXF', '<f8', ()), ('EVENT', 'S4', ()), ('EXISTEVT', 'S4', ()),
             ('NEVNTS', '<i8', ()), ('EVNT_POS', '<i8', ()), ('DELETE', 'S4', ()), ('EXISTDEL', 'S4', ()),
             ('NDELNS', '<i8', ()), ('DEL_POS', '<i8', ()), ('MAP', 'S4', ()), ('EXISTMAP', 'S4', ()),
             ('NMAPS', '<i8', ()), ('MAP_POS', '<i8', ()), ('ACOUST', 'S4', ()), ('EXISTAC', 'S4', ()),
             ('SPATLID', '<i8', ()), ('DISTID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/FATIGUE/TIM2PSD/EVNT', '/NASTRAN/INPUT/FATIGUE/TIM2PSD/DEL',
                 '/NASTRAN/INPUT/FATIGUE/TIM2PSD/MAP']


@register_table
class LOCVAL(object):
    name = 'LOCVAL'
    path = '/NASTRAN/INPUT/LOAD/ACCEL'
    dtype = [('LOC', '<f8', ()), ('VAL', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/ACCEL'
    dtype = [('SID', '<i8', ()), ('CID', '<i8', ()), ('N', '<f8', (3,)), ('DIR', 'S4', ()), ('LOCVAL_LEN', '<i8', ()),
             ('LOCVAL_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/ACCEL/LOCVAL']


@register_table
class GIDL(object):
    name = 'GIDL'
    path = '/NASTRAN/INPUT/LOAD/ACCEL1'
    dtype = [('GIDL', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GTBY(object):
    name = 'GTBY'
    path = '/NASTRAN/INPUT/LOAD/ACCEL1'
    dtype = [('FIRST', '<i8', ()), ('THRU', '<i8', ()), ('BY', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/ACCEL1'
    dtype = [('SID', '<i8', ()), ('CID', '<i8', ()), ('A', '<f8', ()), ('N', '<f8', (3,)), ('GIDL_POS', '<i8', ()),
             ('GIDL_LEN', '<i8', ()), ('GTBY_POS', '<i8', ()), ('GTBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/ACCEL1/GIDL', '/NASTRAN/INPUT/LOAD/ACCEL1/GTBY']


@register_table
class ACSRCE(object):
    name = 'ACSRCE'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('DPHASE', '<i8', ()), ('DELAY', '<i8', ()), ('TC', '<i8', ()),
             ('RHO', '<f8', ()), ('B', '<f8', ()), ('T', '<f8', ()), ('PH', '<f8', ()), ('RC', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEDW(object):
    name = 'AEDW'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('MACH', '<f8', ()), ('SYMXZ', 'S8', ()), ('SYMXY', 'S8', ()), ('UXID', '<i8', ()), ('DMIJ', 'S8', ()),
             ('DMIJI', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEFORCE(object):
    name = 'AEFORCE'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('MACH', '<f8', ()), ('SYMXZ', 'S8', ()), ('SYMXY', 'S8', ()), ('UXID', '<i8', ()), ('MESH', 'S8', ()),
             ('LSET', '<i8', ()), ('DMIK', 'S8', ()), ('PERQ', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEPRESS(object):
    name = 'AEPRESS'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('MACH', '<f8', ()), ('SYMXZ', 'S8', ()), ('SYMXY', 'S8', ()), ('UXID', '<i8', ()), ('DMIJ', 'S8', ()),
             ('DMIJI', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONV(object):
    name = 'CONV'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('EID', '<i8', ()), ('PCONID', '<i8', ()), ('FLMND', '<i8', ()), ('CNTRLND', '<i8', ()),
             ('TA', '<i8', (8,)), ('WT', '<f8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONVM(object):
    name = 'CONVM'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('EID', '<i8', ()), ('PCONID', '<i8', ()), ('FLMND', '<i8', ()), ('CNTMDOT', '<i8', ()),
             ('TA', '<i8', (2,)), ('MDOT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DAREA(object):
    name = 'DAREA'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<i8', ()), ('C', '<i8', ()), ('A', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DELAY(object):
    name = 'DELAY'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<i8', ()), ('C', '<i8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SL(object):
    name = 'SL'
    path = '/NASTRAN/INPUT/LOAD/DLOAD'
    dtype = [('SI', '<f8', ()), ('LI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/DLOAD'
    dtype = [('SID', '<i8', ()), ('S', '<f8', ()), ('SL_LEN', '<i8', ()), ('SL_POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/DLOAD/SL']


@register_table
class DPHASE(object):
    name = 'DPHASE'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<i8', ()), ('C', '<i8', ()), ('TH', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FBODYLD(object):
    name = 'FBODYLD'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('NAMEL', 'S8', ()), ('FBODYSB', 'S8', ()), ('LABEL', 'S64', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FBODYSB(object):
    name = 'FBODYSB'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('NAMES', 'S8', ()), ('GRIDSET', '<i8', ()), ('ELEMSET', '<i8', ()), ('FLAG', '<i8', ()),
             ('LABEL', 'S64', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCE(object):
    name = 'FORCE'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('CID', '<i8', ()), ('F', '<f8', ()), ('N', '<f8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCE1(object):
    name = 'FORCE1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('F', '<f8', ()), ('GI', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCE2(object):
    name = 'FORCE2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('F', '<f8', ()), ('GI', '<i8', (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCEAX(object):
    name = 'FORCEAX'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('RID', '<i8', ()), ('HID1', '<i8', ()), ('HID2', '<i8', ()), ('S', '<f8', ()),
             ('FR', '<f8', ()), ('FP', '<f8', ()), ('FZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRAV(object):
    name = 'GRAV'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('CID', '<i8', ()), ('A', '<f8', ()), ('N', '<f8', (3,)), ('MB', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GUST(object):
    name = 'GUST'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DLOAD', '<i8', ()), ('WG', '<f8', ()), ('X0', '<f8', ()), ('V', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SFACTORS(object):
    name = 'SFACTORS'
    path = '/NASTRAN/INPUT/LOAD/LOAD'
    dtype = [('SI', '<f8', ()), ('LI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/LOAD'
    dtype = [('SID', '<i8', ()), ('S', '<f8', ()), ('SFACTORS_POS', '<i8', ()), ('SFACTORS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/LOAD/SFACTORS']


@register_table
class LOADS(object):
    name = 'LOADS'
    path = '/NASTRAN/INPUT/LOAD/LOADCLID'
    dtype = [('SI', '<f8', ()), ('LIDI', '<i8', ()), ('PARTNAMI', 'S64', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/LOADCLID'
    dtype = [('SID', '<i8', ()), ('S', '<f8', ()), ('LOADS_POS', '<i8', ()), ('LOADS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/LOADCLID/LOADS']


@register_table
class LOADS(object):
    name = 'LOADS'
    path = '/NASTRAN/INPUT/LOAD/LOADCNAM'
    dtype = [('SI', '<f8', ()), ('LOADNAMI', 'S64', ()), ('PARTNAMI', 'S64', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/LOADCNAM'
    dtype = [('SID', '<i8', ()), ('S', '<f8', ()), ('LOADS_POS', '<i8', ()), ('LOADS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/LOADCNAM/LOADS']


@register_table
class LOADS(object):
    name = 'LOADS'
    path = '/NASTRAN/INPUT/LOAD/LOADCSUB'
    dtype = [('SI', '<f8', ()), ('SIDI', '<i8', ()), ('PARTNAMI', 'S64', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/LOADCSUB'
    dtype = [('SID', '<i8', ()), ('S', '<f8', ()), ('LOADS_POS', '<i8', ()), ('LOADS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/LOADCSUB/LOADS']


@register_table
class LSEQ(object):
    name = 'LSEQ'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('LID', '<i8', ()), ('TID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOMAX(object):
    name = 'MOMAX'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('RID', '<i8', ()), ('HID1', '<i8', ()), ('HID2', '<i8', ()), ('S', '<f8', ()),
             ('MR', '<f8', ()), ('MP', '<f8', ()), ('MZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOMENT(object):
    name = 'MOMENT'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('CID', '<i8', ()), ('M', '<f8', ()), ('N', '<f8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOMENT1(object):
    name = 'MOMENT1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('M', '<f8', ()), ('GI', '<i8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOMENT2(object):
    name = 'MOMENT2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('M', '<f8', ()), ('GI', '<i8', (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NLRGAP(object):
    name = 'NLRGAP'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()), ('PLANE', '<i8', ()), ('TABK', '<i8', ()),
             ('TABG', '<i8', ()), ('TABU', '<i8', ()), ('RADIUS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NOLIN1(object):
    name = 'NOLIN1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('S', '<f8', ()), ('GJ', '<i8', ()),
             ('CJ', '<i8', ()), ('T', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NOLIN2(object):
    name = 'NOLIN2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('S', '<f8', ()), ('GJ', '<i8', ()),
             ('CJ', '<i8', ()), ('GK', '<i8', ()), ('CK', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NOLIN3(object):
    name = 'NOLIN3'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('S', '<f8', ()), ('GJ', '<i8', ()),
             ('CJ', '<i8', ()), ('A', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NOLIN4(object):
    name = 'NOLIN4'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('S', '<f8', ()), ('GJ', '<i8', ()),
             ('CJ', '<i8', ()), ('A', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOAD(object):
    name = 'PLOAD'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<f8', ()), ('G', '<i8', (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOAD1(object):
    name = 'PLOAD1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('TYPE', '<i8', ()), ('SCALE', '<i8', ()), ('X1', '<f8', ()),
             ('P1', '<f8', ()), ('X2', '<f8', ()), ('P2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOAD2(object):
    name = 'PLOAD2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<f8', ()), ('EID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOAD4(object):
    name = 'PLOAD4'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('P', '<f8', (4,)), ('G1', '<i8', ()), ('G34', '<i8', ()),
             ('CID', '<i8', ()), ('N', '<f8', (3,)), ('SORL', 'S8', ()), ('LDIR', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOADB3(object):
    name = 'PLOADB3'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('CID', '<i8', ()), ('N', '<f8', (3,)), ('TYPE', '<i8', ()),
             ('SCALE', '<f8', ()), ('P', '<f8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLOADX1(object):
    name = 'PLOADX1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('PA', '<f8', ()), ('PB', '<f8', ()), ('G', '<i8', (2,)),
             ('THETA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRESAX(object):
    name = 'PRESAX'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('P', '<f8', ()), ('RID', '<i8', (2,)), ('PHI', '<f8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QBDY1(object):
    name = 'QBDY1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('Q0', '<f8', ()), ('EID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QBDY2(object):
    name = 'QBDY2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('EID', '<i8', ()), ('Q0', '<f8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QBDY3(object):
    name = 'QBDY3'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('Q0', '<f8', ()), ('CNTRLND', '<i8', ()), ('EID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QHBDY(object):
    name = 'QHBDY'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('FLAG', '<i8', ()), ('Q0', '<f8', ()), ('AF', '<f8', ()), ('G', '<i8', (8,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QVECT(object):
    name = 'QVECT'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('Q0', '<f8', ()), ('TSOUR', '<f8', ()), ('CE', '<i8', ()), ('TID', '<i8', (3,)),
             ('E', '<f8', (3,)), ('CNTRLND', '<i8', ()), ('EID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QVOL(object):
    name = 'QVOL'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('QVOL', '<f8', ()), ('CNTRLND', '<i8', ()), ('EID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RADBC(object):
    name = 'RADBC'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('EID', '<i8', ()), ('FAMB', '<f8', ()), ('CNTRLND', '<i8', ()), ('NODAMB', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RFORCE(object):
    name = 'RFORCE'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('CID', '<i8', ()), ('A', '<f8', ()), ('R', '<f8', (3,)),
             ('METHOD', '<i8', ()), ('RACC', '<f8', ()), ('MB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RLOAD1(object):
    name = 'RLOAD1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('DELAY', '<i8', ()), ('DPHASE', '<i8', ()), ('TC', '<i8', ()),
             ('TD', '<i8', ()), ('TYPE', '<i8', ()), ('T', '<f8', ()), ('PH', '<f8', ()), ('RC', '<f8', ()),
             ('RD', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RLOAD2(object):
    name = 'RLOAD2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('DELAY', '<i8', ()), ('DPHASE', '<i8', ()), ('TB', '<i8', ()),
             ('TP', '<i8', ()), ('TYPE', '<i8', ()), ('T', '<f8', ()), ('PH', '<f8', ()), ('RB', '<f8', ()),
             ('RP', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOAD(object):
    name = 'SLOAD'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('F', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOADN1(object):
    name = 'SLOADN1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('D', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SIGMA(object):
    name = 'SIGMA'
    path = '/NASTRAN/INPUT/LOAD/TF'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('A0I', '<f8', ()), ('A1I', '<f8', ()), ('A2I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/LOAD/TF'
    dtype = [('SID', '<i8', ()), ('GD', '<i8', ()), ('CD', '<i8', ()), ('B0', '<f8', ()), ('B1', '<f8', ()),
             ('B2', '<f8', ()), ('SIGMA_POS', '<i8', ()), ('SIGMA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/LOAD/TF/SIGMA']


@register_table
class TIC(object):
    name = 'TIC'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('G', '<i8', ()), ('C', '<i8', ()), ('U0', '<f8', ()), ('V0', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TLOAD1(object):
    name = 'TLOAD1'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('DELAY', '<i8', ()), ('TYPE', '<i8', ()), ('TID', '<i8', ()),
             ('U0', '<f8', ()), ('V0', '<f8', ()), ('T', '<f8', ()), ('F', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TLOAD2(object):
    name = 'TLOAD2'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('DAREA', '<i8', ()), ('DELAY', '<i8', ()), ('TYPE', '<i8', ()), ('T1', '<f8', ()),
             ('T2', '<f8', ()), ('F', '<f8', ()), ('P', '<f8', ()), ('C', '<f8', ()), ('B', '<f8', ()),
             ('U0', '<f8', ()), ('V0', '<f8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TTEMP(object):
    name = 'TTEMP'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('GROUPID', '<i8', ()), ('TID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class UNBALNC(object):
    name = 'UNBALNC'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('SID', '<i8', ()), ('MTABLE', '<i8', ()), ('MASS', '<f8', ()), ('GRID', '<i8', ()), ('X1', '<f8', ()),
             ('X2', '<f8', ()), ('X3', '<f8', ()), ('ROFFTAB', '<i8', ()), ('ROFFSET', '<f8', ()), ('THETA', '<f8', ()),
             ('ZOFFTAB', '<i8', ()), ('ZOFFSET', '<f8', ()), ('TON', '<f8', ()), ('TOFF', '<f8', ()),
             ('CFLAG', 'S8', ()), ('UFT1', '<i8', ()), ('UFT2', '<i8', ()), ('UFT3', '<i8', ()), ('UFR1', '<i8', ()),
             ('UFR2', '<i8', ()), ('UFR3', '<i8', ()), ('MCT1', '<i8', ()), ('MCT2', '<i8', ()), ('MCT3', '<i8', ()),
             ('MCR1', '<i8', ()), ('MCR2', '<i8', ()), ('MCR3', '<i8', ()), ('SCR1', '<i8', ()), ('SCR2', '<i8', ()),
             ('SCR3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WETLOAD(object):
    name = 'WETLOAD'
    path = '/NASTRAN/INPUT/LOAD'
    dtype = [('WLID', '<i8', ()), ('WSID', '<i8', ()), ('SERV_ID', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CREEP(object):
    name = 'CREEP'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('T0', '<f8', ()), ('EXP', '<f8', ()), ('FORM', '<i8', ()), ('TIDKP', '<i8', ()),
             ('TIDCP', '<i8', ()), ('TIDCS', '<i8', ()), ('THRESH', '<f8', ()), ('TYPE', '<i8', ()),
             ('AG', '<f8', (7,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT1(object):
    name = 'MAT1'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('E', '<f8', ()), ('G', '<f8', ()), ('NU', '<f8', ()), ('RHO', '<f8', ()),
             ('A', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('ST', '<f8', ()), ('SC', '<f8', ()),
             ('SS', '<f8', ()), ('MCSID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT10(object):
    name = 'MAT10'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('BULK', '<f8', ()), ('RHO', '<f8', ()), ('C', '<f8', ()), ('GE', '<f8', ()),
             ('POROC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT2(object):
    name = 'MAT2'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('G11', '<f8', ()), ('G12', '<f8', ()), ('G13', '<f8', ()), ('G22', '<f8', ()),
             ('G23', '<f8', ()), ('G33', '<f8', ()), ('RHO', '<f8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()),
             ('A12', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('ST', '<f8', ()), ('SC', '<f8', ()),
             ('SS', '<f8', ()), ('MCSID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT3(object):
    name = 'MAT3'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('EX', '<f8', ()), ('ETH', '<f8', ()), ('EZ', '<f8', ()), ('NUXTH', '<f8', ()),
             ('NUTHZ', '<f8', ()), ('NUZX', '<f8', ()), ('RHO', '<f8', ()), ('GZX', '<f8', ()), ('AX', '<f8', ()),
             ('ATH', '<f8', ()), ('AZ', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT4(object):
    name = 'MAT4'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('K', '<f8', ()), ('CP', '<f8', ()), ('RHO', '<f8', ()), ('H', '<f8', ()),
             ('MU', '<f8', ()), ('HGEN', '<f8', ()), ('REFENTH', '<f8', ()), ('TCH', '<f8', ()), ('TDELTA', '<f8', ()),
             ('QLAT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT5(object):
    name = 'MAT5'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('K', '<f8', (6,)), ('CP', '<f8', ()), ('RHO', '<f8', ()), ('HGEN', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT8(object):
    name = 'MAT8'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('E1', '<f8', ()), ('E2', '<f8', ()), ('NU12', '<f8', ()), ('G12', '<f8', ()),
             ('G1Z', '<f8', ()), ('G2Z', '<f8', ()), ('RHO', '<f8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()),
             ('TREF', '<f8', ()), ('XT', '<f8', ()), ('XC', '<f8', ()), ('YT', '<f8', ()), ('YC', '<f8', ()),
             ('S', '<f8', ()), ('GE', '<f8', ()), ('F12', '<f8', ()), ('STRN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MAT9(object):
    name = 'MAT9'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('G', '<f8', (21,)), ('RHO', '<f8', ()), ('A', '<f8', (6,)), ('TREF', '<f8', ()),
             ('GE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATDIGI(object):
    name = 'MATDIGI'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('UDID', '<i8', ()), ('K', '<f8', ()), ('NU', '<f8', ()), ('E', '<f8', ()),
             ('G', '<f8', ()), ('RHO', '<f8', ()), ('IORTHO', '<i8', ()), ('IFAIL', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ANISO(object):
    name = 'ANISO'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('R11', '<f8', ()),
             ('R22', '<f8', ()), ('R33', '<f8', ()), ('R12', '<f8', ()), ('R23', '<f8', ()), ('R31', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHABOCHE(object):
    name = 'CHABOCHE'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('R0', '<f8', ()),
             ('RINF', '<f8', ()), ('B', '<f8', ()), ('C', '<f8', ()), ('GAM', '<f8', ()), ('KAP', '<f8', ()),
             ('N', '<i8', ()), ('QM', '<f8', ()), ('MU', '<f8', ()), ('ETA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GURSON(object):
    name = 'GURSON'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('Q1', '<f8', ()),
             ('Q2', '<f8', ()), ('INITIAL', '<f8', ()), ('CRITICAL', '<f8', ()), ('FAILURE', '<f8', ()),
             ('NUCL', 'S8', ()), ('MEAN', '<f8', ()), ('SDEV', '<f8', ()), ('NFRAC', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IMPCREEP(object):
    name = 'IMPCREEP'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('VMISES', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class JHCOOK(object):
    name = 'JHCOOK'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('A', '<f8', ()),
             ('B', '<f8', ()), ('N', '<f8', ()), ('C', '<f8', ()), ('M', '<f8', ()), ('TMELT', '<f8', ()),
             ('TROOM', '<f8', ()), ('E0DOT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class KUMAR(object):
    name = 'KUMAR'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('B0', '<f8', ()),
             ('A', '<f8', ()), ('B1', '<f8', ()), ('B2', '<f8', ()), ('B3', '<f8', ()), ('TB0', '<i8', ()),
             ('TA', '<i8', ()), ('N', '<f8', ()), ('B4', '<f8', ()), ('B5', '<f8', ()), ('B6', '<f8', ()),
             ('TN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ORNL(object):
    name = 'ORNL'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('OPTION', 'S8', ()),
             ('YC10', '<f8', ()), ('TID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRESS(object):
    name = 'PRESS'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('OPTION', 'S8', ()),
             ('ALPHA', '<f8', ()), ('BETA', '<f8', ()), ('CRACKS', '<f8', ()), ('SOFTEN', '<f8', ()),
             ('CRUSHS', '<f8', ()), ('SRFAC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIMARY(object):
    name = 'PRIMARY'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PWRLAW(object):
    name = 'PWRLAW'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('A', '<f8', ()),
             ('M', '<f8', ()), ('B', '<f8', ()), ('N', '<f8', ()), ('S0E0', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class REFFECT(object):
    name = 'REFFECT'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('OPTION', 'S8', ()),
             ('RTID', '<i8', ()), ('C', '<f8', ()), ('P', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class UNITS(object):
    name = 'UNITS'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('UOPT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VPARAM(object):
    name = 'VPARAM'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('NVP', '<i8', ()),
             ('VP1', '<f8', ()), ('VP2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class YLDOPT(object):
    name = 'YLDOPT'
    path = '/NASTRAN/INPUT/MATERIAL/MATEP'
    dtype = [('MID', '<i8', ()), ('FORM', 'S8', ()), ('Y0', '<f8', ()), ('FID', '<i8', ()), ('RYIELD', 'S8', ()),
             ('WKHARD', 'S8', ()), ('METHOD', 'S8', ()), ('H', '<f8', ()), ('ITYPE', '<i8', ()), ('OPTION', 'S8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATF_DATA(object):
    name = 'MATF_DATA'
    path = '/NASTRAN/INPUT/MATERIAL/MATF'
    dtype = [('CRI', 'S4', ()), ('CRITERIA', '<i8', ()), ('XT', '<f8', ()), ('XC', '<f8', ()), ('YT', '<f8', ()),
             ('YC', '<f8', ()), ('ZT', '<f8', ()), ('ZC', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
             ('SZX', '<f8', ()), ('FIND', '<f8', ()), ('FXY', '<f8', ()), ('FYZ', '<f8', ()), ('FZX', '<f8', ()),
             ('EXT', '<f8', ()), ('EXC', '<f8', ()), ('EYT', '<f8', ()), ('EYC', '<f8', ()), ('EZT', '<f8', ()),
             ('EZC', '<f8', ()), ('GXY', '<f8', ()), ('GYZ', '<f8', ()), ('GZX', '<f8', ()), ('PF', 'S4', ()),
             ('A1', '<f8', ()), ('A2', '<f8', ()), ('A3', '<f8', ()), ('A4', '<f8', ()), ('A5', '<f8', ()),
             ('IC1', '<i8', ()), ('IC2', '<i8', ()), ('IC3', '<i8', ()), ('IC4', '<i8', ()), ('IC5', '<i8', ()),
             ('IC6', '<i8', ()), ('IC7', '<i8', ()), ('IC8', '<i8', ()), ('IC9', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIMARY(object):
    name = 'PRIMARY'
    path = '/NASTRAN/INPUT/MATERIAL/MATF'
    dtype = [('MID', '<i8', ()), ('ITYPE', '<i8', ()), ('SB', '<f8', ()), ('NCRI', '<i8', ()), ('D_POS', '<i8', ()),
             ('D_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATF1(object):
    name = 'MATF1'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE', '<i8', ()), ('TG', '<i8', ()), ('TNU', '<i8', ()), ('TRHO', '<i8', ()),
             ('TA', '<i8', ()), ('TGE', '<i8', ()), ('TST', '<i8', ()), ('TSC', '<i8', ()), ('TSS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATG(object):
    name = 'MATG'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('IDMEM', '<i8', ()), ('BEHAV', '<i8', ()), ('TABLD', '<i8', ()),
             ('TABLU1', '<i8', ()), ('TABLU2', '<i8', ()), ('TABLU3', '<i8', ()), ('TABLU4', '<i8', ()),
             ('TABLU5', '<i8', ()), ('TABLU6', '<i8', ()), ('TABLU7', '<i8', ()), ('TABLU8', '<i8', ()),
             ('TABLU9', '<i8', ()), ('TABLU10', '<i8', ()), ('YPRS', '<f8', ()), ('EPL', '<f8', ()), ('GPL', '<f8', ()),
             ('GAP', '<f8', ()), ('TABYPRS', '<i8', ()), ('TABEPL', '<i8', ()), ('TABGPL', '<i8', ()),
             ('TABGAP', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ABOYCE_GENT(object):
    name = 'ABOYCE_GENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('NOT', '<i8', ()), ('K', '<f8', ()), ('RHO', '<f8', ()),
             ('TEXP', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('NKT', '<f8', ()), ('NE', '<f8', ()),
             ('IM', '<f8', ()), ('DI', '<f8', (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GHEM(object):
    name = 'GHEM'
    path = '/NASTRAN/INPUT/MATERIAL/MATHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('NOT', '<i8', ()), ('K', '<f8', ()), ('RHO', '<f8', ()),
             ('TEXP', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOONEY(object):
    name = 'MOONEY'
    path = '/NASTRAN/INPUT/MATERIAL/MATHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('NOT', '<i8', ()), ('K', '<f8', ()), ('RHO', '<f8', ()),
             ('TEXP', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('C10', '<f8', ()), ('C01', '<f8', ()),
             ('D1', '<f8', ()), ('TAB1', '<i8', ()), ('TAB2', '<i8', ()), ('TAB3', '<i8', ()), ('TAB4', '<i8', ()),
             ('TABD', '<i8', ()), ('C20', '<f8', ()), ('C11', '<f8', ()), ('C02', '<f8', ()), ('D2', '<f8', ()),
             ('NA', '<i8', ()), ('ND', '<i8', ()), ('C30', '<f8', ()), ('C21', '<f8', ()), ('C12', '<f8', ()),
             ('C03', '<f8', ()), ('D3', '<f8', ()), ('C40', '<f8', ()), ('C31', '<f8', ()), ('C22', '<f8', ()),
             ('C13', '<f8', ()), ('C04', '<f8', ()), ('D4', '<f8', ()), ('C50', '<f8', ()), ('C41', '<f8', ()),
             ('C32', '<f8', ()), ('C23', '<f8', ()), ('C14', '<f8', ()), ('C05', '<f8', ()), ('D5', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class OGDEN_FOAM(object):
    name = 'OGDEN_FOAM'
    path = '/NASTRAN/INPUT/MATERIAL/MATHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('NOT', '<i8', ()), ('K', '<f8', ()), ('RHO', '<f8', ()),
             ('TEXP', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('MU1', '<f8', ()), ('ALPHA1', '<f8', ()),
             ('BETA1', '<f8', ()), ('MU2', '<f8', ()), ('ALPHA2', '<f8', ()), ('BETA2', '<f8', ()), ('MU3', '<f8', ()),
             ('ALPHA3', '<f8', ()), ('BETA3', '<f8', ()), ('MU4', '<f8', ()), ('ALPHA4', '<f8', ()),
             ('BETA4', '<f8', ()), ('MU5', '<f8', ()), ('ALPHA5', '<f8', ()), ('BETA5', '<f8', ()), ('DI', '<f8', (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATHP'
    dtype = [('MID', '<i8', ()), ('A10', '<f8', ()), ('A01', '<f8', ()), ('D1', '<f8', ()), ('RHO', '<f8', ()),
             ('ALPHA', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('SF', '<i8', ()), ('NA', '<i8', ()),
             ('ND', '<i8', ()), ('KP', '<f8', ()), ('A20', '<f8', ()), ('A11', '<f8', ()), ('A02', '<f8', ()),
             ('D2', '<f8', ()), ('A30', '<f8', ()), ('A21', '<f8', ()), ('A12', '<f8', ()), ('A03', '<f8', ()),
             ('D3', '<f8', ()), ('A40', '<f8', ()), ('A31', '<f8', ()), ('A22', '<f8', ()), ('A13', '<f8', ()),
             ('A04', '<f8', ()), ('D4', '<f8', ()), ('A50', '<f8', ()), ('A41', '<f8', ()), ('A32', '<f8', ()),
             ('A23', '<f8', ()), ('A14', '<f8', ()), ('A05', '<f8', ()), ('D5', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TABLE(object):
    name = 'TABLE'
    path = '/NASTRAN/INPUT/MATERIAL/MATHP'
    dtype = [('TAB1', '<i8', ()), ('TAB2', '<i8', ()), ('TAB3', '<i8', ()), ('TAB4', '<i8', ()), ('TAB5', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATORT(object):
    name = 'MATORT'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('E1', '<f8', ()), ('E2', '<f8', ()), ('E3', '<f8', ()), ('NU12', '<f8', ()),
             ('NU23', '<f8', ()), ('NU31', '<f8', ()), ('RHO', '<f8', ()), ('G12', '<f8', ()), ('G23', '<f8', ()),
             ('G31', '<f8', ()), ('A1', '<f8', ()), ('A2', '<f8', ()), ('A3', '<f8', ()), ('TREF', '<f8', ()),
             ('GE', '<f8', ()), ('IYLD', '<i8', ()), ('IHARD', '<i8', ()), ('SY', '<f8', ()), ('SORNL', '<f8', ()),
             ('Y1', '<f8', ()), ('Y2', '<f8', ()), ('Y3', '<f8', ()), ('YSHR1', '<f8', ()), ('YSHR2', '<f8', ()),
             ('YSHR3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATPE1(object):
    name = 'MATPE1'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('KEY', '<i8', ()), ('MAT1', '<i8', ()), ('MAT10', '<i8', ()), ('BIOT', '<f8', ()),
             ('VISC', '<f8', ()), ('GAMMA', '<f8', ()), ('PRANDTL', '<f8', ()), ('POR', '<f8', ()), ('TOR', '<f8', ()),
             ('AFR', '<f8', ()), ('VLE', '<f8', ()), ('TLE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATS1(object):
    name = 'MATS1'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TID', '<i8', ()), ('TYPE', '<i8', ()), ('H', '<f8', ()), ('YF', '<i8', ()),
             ('HR', '<i8', ()), ('LIMIT1', '<f8', ()), ('LIMIT2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATS3(object):
    name = 'MATS3'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TEX', '<i8', ()), ('TETH', '<i8', ()), ('TEZ', '<i8', ()), ('TNUXTH', '<i8', ()),
             ('TNUTHZ', '<i8', ()), ('TNUZX', '<i8', ()), ('TRHO', '<i8', ()), ('TGZX', '<i8', ()), ('TAX', '<i8', ()),
             ('TATH', '<i8', ()), ('TAZ', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATS8(object):
    name = 'MATS8'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE1', '<i8', ()), ('TE2', '<i8', ()), ('TNU12', '<i8', ()), ('TG12', '<i8', ()),
             ('TG1Z', '<i8', ()), ('TG2Z', '<i8', ()), ('TRHO', '<i8', ()), ('TGZX', '<i8', ()), ('TA1', '<i8', ()),
             ('TA2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATSMA(object):
    name = 'MATSMA'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('T0', '<f8', ()), ('EL', '<f8', ()), ('AUSPROP', '<f8', (8,)),
             ('MARPROP', '<f8', (8,)), ('THMOMECH', '<f8', (12,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATSORT(object):
    name = 'MATSORT'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE1', '<i8', ()), ('TE2', '<i8', ()), ('TE3', '<i8', ()), ('TNU12', '<i8', ()),
             ('TNU23', '<i8', ()), ('TNU31', '<i8', ()), ('TRHO', '<i8', ()), ('TG12', '<i8', ()), ('TG23', '<i8', ()),
             ('TG31', '<i8', ()), ('TA1', '<i8', ()), ('TA2', '<i8', ()), ('TA3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT1(object):
    name = 'MATT1'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE', '<i8', ()), ('TG', '<i8', ()), ('TNU', '<i8', ()), ('TRHO', '<i8', ()),
             ('TA', '<i8', ()), ('TGE', '<i8', ()), ('TST', '<i8', ()), ('TSC', '<i8', ()), ('TSS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT2(object):
    name = 'MATT2'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TG11', '<i8', ()), ('TG12', '<i8', ()), ('TG13', '<i8', ()), ('TG22', '<i8', ()),
             ('TG23', '<i8', ()), ('TG33', '<i8', ()), ('TRHO', '<i8', ()), ('TA1', '<i8', ()), ('TA2', '<i8', ()),
             ('TA3', '<i8', ()), ('TGE', '<i8', ()), ('TST', '<i8', ()), ('TSC', '<i8', ()), ('TSS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT3(object):
    name = 'MATT3'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TID', '<i8', (15,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT4(object):
    name = 'MATT4'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TK', '<i8', ()), ('TCP', '<i8', ()), ('TH', '<i8', ()), ('TMU', '<i8', ()),
             ('THGEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT5(object):
    name = 'MATT5'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TK', '<i8', (6,)), ('TCP', '<i8', ()), ('THGEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT8(object):
    name = 'MATT8'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE1', '<i8', ()), ('TE2', '<i8', ()), ('TNU12', '<i8', ()), ('TG12', '<i8', ()),
             ('TG1Z', '<i8', ()), ('TG2Z', '<i8', ()), ('TRHO', '<i8', ()), ('TA1', '<i8', ()), ('TA2', '<i8', ()),
             ('TXT', '<i8', ()), ('TXC', '<i8', ()), ('TYT', '<i8', ()), ('TYC', '<i8', ()), ('TS', '<i8', ()),
             ('TGE', '<i8', ()), ('TF12', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATT9(object):
    name = 'MATT9'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TG', '<i8', (21,)), ('TRHO', '<i8', ()), ('TA', '<i8', (6,)), ('TGE', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHABOCHE(object):
    name = 'CHABOCHE'
    path = '/NASTRAN/INPUT/MATERIAL/MATTEP'
    dtype = [('R0', '<i8', ()), ('RINF', '<i8', ()), ('B', '<i8', ()), ('C', '<i8', ()), ('GAM', '<i8', ()),
             ('KAP', '<i8', ()), ('N', '<i8', ()), ('QM', '<i8', ()), ('MU', '<i8', ()), ('ETA', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATTEP'
    dtype = [('MID', '<i8', ()), ('TY0', '<i8', ()), ('TFID', '<i8', ()), ('TH', '<i8', ()), ('TYC10', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class POWER(object):
    name = 'POWER'
    path = '/NASTRAN/INPUT/MATERIAL/MATTEP'
    dtype = [('A', '<i8', ()), ('M', '<i8', ()), ('B', '<i8', ()), ('N', '<i8', ()), ('S0E0', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATTF_DATA(object):
    name = 'MATTF_DATA'
    path = '/NASTRAN/INPUT/MATERIAL/MATTF'
    dtype = [('KIND', '<i8', ()), ('CRITERIA', '<i8', ()), ('TIDXT', '<i8', ()), ('TIDXC', '<i8', ()),
             ('TIDYT', '<i8', ()), ('TIDYC', '<i8', ()), ('TIDZT', '<i8', ()), ('TIDZC', '<i8', ()),
             ('TIDSXY', '<i8', ()), ('TIDSYZ', '<i8', ()), ('TIDSZX', '<i8', ()), ('TIDFIND', '<i8', ()),
             ('TIDFXY', '<i8', ()), ('TIDFYZ', '<i8', ()), ('TIDFZX', '<i8', ()), ('TIDEXT', '<i8', ()),
             ('TIDEXC', '<i8', ()), ('TIDEYT', '<i8', ()), ('TIDEYC', '<i8', ()), ('TIDEZT', '<i8', ()),
             ('TIDEZC', '<i8', ()), ('TIDGXY', '<i8', ()), ('TIDGYZ', '<i8', ()), ('TIDGZX', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRIMARY(object):
    name = 'PRIMARY'
    path = '/NASTRAN/INPUT/MATERIAL/MATTF'
    dtype = [('MID', '<i8', ()), ('NCRI', '<i8', ()), ('D_POS', '<i8', ()), ('D_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATTG(object):
    name = 'MATTG'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('IDYM', '<i8', ()), ('IDVM', '<i8', ()), ('IDDM', '<i8', ()), ('IDLD', '<i8', ()),
             ('IDU', '<i8', (10,)), ('IDYPR', '<i8', ()), ('IDEPL', '<i8', ()), ('IDGPL', '<i8', ()),
             ('IDGAP', '<i8', ()), ('TABEPL', '<i8', ()), ('TABGPL', '<i8', ()), ('TABGAP', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ABOYCE_GENT(object):
    name = 'ABOYCE_GENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATTHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('K', '<i8', ()), ('RHO', '<i8', ()), ('TEXP', '<i8', ()),
             ('GE', '<i8', ()), ('NKT', '<i8', ()), ('NE', '<i8', ()), ('IM', '<i8', ()), ('DI', '<i8', (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATTHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('K', '<i8', ()), ('RHO', '<i8', ()), ('TEXP', '<i8', ()),
             ('GE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOONEY(object):
    name = 'MOONEY'
    path = '/NASTRAN/INPUT/MATERIAL/MATTHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('K', '<i8', ()), ('RHO', '<i8', ()), ('TEXP', '<i8', ()),
             ('GE', '<i8', ()), ('C10', '<i8', ()), ('C01', '<i8', ()), ('D1', '<i8', ()), ('C20', '<i8', ()),
             ('C11', '<i8', ()), ('C02', '<i8', ()), ('D2', '<i8', ()), ('C30', '<i8', ()), ('C21', '<i8', ()),
             ('C12', '<i8', ()), ('C03', '<i8', ()), ('D3', '<i8', ()), ('C40', '<i8', ()), ('C31', '<i8', ()),
             ('C22', '<i8', ()), ('C13', '<i8', ()), ('C04', '<i8', ()), ('D4', '<i8', ()), ('C50', '<i8', ()),
             ('C41', '<i8', ()), ('C32', '<i8', ()), ('C23', '<i8', ()), ('C14', '<i8', ()), ('C05', '<i8', ()),
             ('D5', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class OGDEN_FOAM(object):
    name = 'OGDEN_FOAM'
    path = '/NASTRAN/INPUT/MATERIAL/MATTHE'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('K', '<i8', ()), ('RHO', '<i8', ()), ('TEXP', '<i8', ()),
             ('GE', '<i8', ()), ('MU1', '<i8', ()), ('ALPHA1', '<i8', ()), ('BETA1', '<i8', ()), ('MU2', '<i8', ()),
             ('ALPHA2', '<i8', ()), ('BETA2', '<i8', ()), ('MU3', '<i8', ()), ('ALPHA3', '<i8', ()),
             ('BETA3', '<i8', ()), ('MU4', '<i8', ()), ('ALPHA4', '<i8', ()), ('BETA4', '<i8', ()), ('MU5', '<i8', ()),
             ('ALPHA5', '<i8', ()), ('BETA5', '<i8', ()), ('DI', '<i8', (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATTORT(object):
    name = 'MATTORT'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TE1', '<i8', ()), ('TE2', '<i8', ()), ('TE3', '<i8', ()), ('TNU12', '<i8', ()),
             ('TNU23', '<i8', ()), ('TNU31', '<i8', ()), ('TG12', '<i8', ()), ('TG23', '<i8', ()), ('TG31', '<i8', ()),
             ('TA1', '<i8', ()), ('TA2', '<i8', ()), ('TA3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATTUSR(object):
    name = 'MATTUSR'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('TRHO', '<i8', ()), ('TA1', '<i8', ()), ('TA2', '<i8', ()), ('TA3', '<i8', ()),
             ('TGE', '<i8', ()), ('TST', '<i8', ()), ('TSC', '<i8', ()), ('TSS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NARA(object):
    name = 'NARA'
    path = '/NASTRAN/INPUT/MATERIAL/MATTVE'
    dtype = [('WI', '<f8', ()), ('TI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class POWER(object):
    name = 'POWER'
    path = '/NASTRAN/INPUT/MATERIAL/MATTVE'
    dtype = [('CI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class WLF(object):
    name = 'WLF'
    path = '/NASTRAN/INPUT/MATERIAL/MATTVE'
    dtype = [('A1', '<f8', ()), ('A2', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATERIAL/MATTVE'
    dtype = [('MID', '<i8', ()), ('FUNCT', 'S4', ()), ('IFUNCT', '<i8', ()), ('RT', '<f8', ()), ('ENER', '<f8', ()),
             ('FRACT', '<f8', ()), ('TDIF', '<f8', ()), ('TREF', '<f8', ()), ('NP', '<i8', ()), ('POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATERIAL/MATTVE/NARA', '/NASTRAN/INPUT/MATERIAL/MATTVE/POWER',
                 '/NASTRAN/INPUT/MATERIAL/MATTVE/WLF']


@register_table
class CDATA(object):
    name = 'CDATA'
    path = '/NASTRAN/INPUT/MATERIAL/MATUDS'
    dtype = [('CDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDATA(object):
    name = 'IDATA'
    path = '/NASTRAN/INPUT/MATERIAL/MATUDS'
    dtype = [('IDATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class RDATA(object):
    name = 'RDATA'
    path = '/NASTRAN/INPUT/MATERIAL/MATUDS'
    dtype = [('RDATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATERIAL/MATUDS'
    dtype = [('NWD', '<i8', ()), ('MID', '<i8', ()), ('MTYPE', 'S8', ()), ('GROUP', 'S8', ()), ('UNAME', 'S8', ()),
             ('IOPT', '<i8', ()), ('NINTG', '<i8', ()), ('NREAL', '<i8', ()), ('NCHAR', '<i8', ()), ('IPOS', '<i8', ()),
             ('RPOS', '<i8', ()), ('CPOS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATERIAL/MATUDS/CDATA', '/NASTRAN/INPUT/MATERIAL/MATUDS/IDATA',
                 '/NASTRAN/INPUT/MATERIAL/MATUDS/RDATA']


@register_table
class MATUSR(object):
    name = 'MATUSR'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('IPREF', '<i8', ()), ('IKINEM', '<i8', ()), ('RHO', '<f8', ()), ('A1', '<f8', ()),
             ('A2', '<f8', ()), ('A12', '<f8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('ST', '<f8', ()),
             ('SC', '<f8', ()), ('SS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATVE_IDENT(object):
    name = 'MATVE_IDENT'
    path = '/NASTRAN/INPUT/MATERIAL/MATVE'
    dtype = [('MID', '<i8', ()), ('MODEL', 'S4', ()), ('IMODEL', '<i8', ()), ('POS', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATVE_NONORTHO(object):
    name = 'MATVE_NONORTHO'
    path = '/NASTRAN/INPUT/MATERIAL/MATVE'
    dtype = [('ALPHAS', '<f8', ()), ('ALPHA1', '<f8', ()), ('WD1', '<f8', ()), ('TD1', '<f8', ()), ('WV1', '<f8', ()),
             ('TV1', '<f8', ()), ('WD2', '<f8', ()), ('TD2', '<f8', ()), ('WD3', '<f8', ()), ('TD3', '<f8', ()),
             ('WD4', '<f8', ()), ('TD4', '<f8', ()), ('WD5', '<f8', ()), ('TD5', '<f8', ()), ('WV2', '<f8', ()),
             ('TV2', '<f8', ()), ('WV3', '<f8', ()), ('TV3', '<f8', ()), ('WV4', '<f8', ()), ('TV4', '<f8', ()),
             ('WV5', '<f8', ()), ('TV5', '<f8', ()), ('NUMDEV', '<i8', ()), ('NUMVOL', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MATVE_ORTHO(object):
    name = 'MATVE_ORTHO'
    path = '/NASTRAN/INPUT/MATERIAL/MATVE'
    dtype = [('TD1', '<f8', ()), ('EXX1', '<f8', ()), ('EYY1', '<f8', ()), ('EZZ1', '<f8', ()), ('VXY1', '<f8', ()),
             ('VYZ1', '<f8', ()), ('VZX1', '<f8', ()), ('GXY1', '<f8', ()), ('GYZ1', '<f8', ()), ('GZX1', '<f8', ()),
             ('TD2', '<f8', ()), ('EXX2', '<f8', ()), ('EYY2', '<f8', ()), ('EZZ2', '<f8', ()), ('VXY2', '<f8', ()),
             ('VYZ2', '<f8', ()), ('VZX2', '<f8', ()), ('GXY2', '<f8', ()), ('GYZ2', '<f8', ()), ('GZX2', '<f8', ()),
             ('TD3', '<f8', ()), ('EXX3', '<f8', ()), ('EYY3', '<f8', ()), ('EZZ3', '<f8', ()), ('VXY3', '<f8', ()),
             ('VYZ3', '<f8', ()), ('VZX3', '<f8', ()), ('GXY3', '<f8', ()), ('GYZ3', '<f8', ()), ('GZX3', '<f8', ()),
             ('TD4', '<f8', ()), ('EXX4', '<f8', ()), ('EYY4', '<f8', ()), ('EZZ4', '<f8', ()), ('VXY4', '<f8', ()),
             ('VYZ4', '<f8', ()), ('VZX4', '<f8', ()), ('GXY4', '<f8', ()), ('GYZ4', '<f8', ()), ('GZX4', '<f8', ()),
             ('TD5', '<f8', ()), ('EXX5', '<f8', ()), ('EYY5', '<f8', ()), ('EZZ5', '<f8', ()), ('VXY5', '<f8', ()),
             ('VYZ5', '<f8', ()), ('VZX5', '<f8', ()), ('GXY5', '<f8', ()), ('GYZ5', '<f8', ()), ('GZX5', '<f8', ()),
             ('NUMDEV', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GENERAL(object):
    name = 'GENERAL'
    path = '/NASTRAN/INPUT/MATERIAL/MATVP'
    dtype = [('MID', '<i8', ()), ('A', '<f8', ()), ('IT3D', '<i8', ()), ('M', '<f8', ()), ('N', '<f8', ()),
             ('P', '<f8', ()), ('Q', '<f8', ()), ('IUSER', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ANAND(object):
    name = 'ANAND'
    path = '/NASTRAN/INPUT/MATERIAL/MATVP'
    dtype = [('MID', '<i8', ()), ('ANAND', 'S8', ()), ('PREXF', '<f8', ()), ('ACTEN', '<f8', ()), ('MULST', '<f8', ()),
             ('STNRT', '<f8', ()), ('SATCO', '<f8', ()), ('STNSA', '<f8', ()), ('HRCN', '<f8', ()),
             ('STNHR', '<f8', ()), ('DEFRS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MCOHE(object):
    name = 'MCOHE'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('MID', '<i8', ()), ('MODEL', '<i8', ()), ('DEACT', '<i8', ()), ('TID', '<i8', ()), ('COHE', '<f8', ()),
             ('CRTOD', '<f8', ()), ('MAXOD', '<f8', ()), ('SNSR', '<f8', ()), ('EXP', '<f8', ()), ('VED', '<f8', ()),
             ('RRRD', '<f8', ()), ('SFC', '<f8', ()), ('SNER', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NLMOPTS(object):
    name = 'NLMOPTS'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('CREEP', 'S4', ()), ('VALC1', '<i8', ()), ('VALC2', '<i8', ()), ('VALC3', '<i8', ()),
             ('VALC4', '<i8', ()), ('ASSM', 'S4', ()), ('VALA', 'S8', ()), ('TSHEAR', 'S4', ()), ('VALT', 'S8', ()),
             ('LRGSTRN', 'S4', ()), ('VALL', '<i8', ()), ('HEMI', 'S4', ()), ('VIEW', '<i8', ()), ('NPIXEL', '<i8', ()),
             ('NDIV', '<i8', ()), ('CUTOFF', '<f8', ()), ('FRACTION', '<f8', ()), ('FACCNT', '<i8', ()),
             ('FACTOL', '<f8', ()), ('SPROPMAP', 'S4', ()), ('PROPMAP', '<i8', ()), ('PROPBEH', 'S4', ()),
             ('DIRECT', '<i8', ()), ('THICKOP', '<f8', ()), ('IPRINT', '<i8', ()), ('TEMPP', 'S4', ()),
             ('TMPDIST', 'S4', ()), ('TEMGO', 'S4', ()), ('VMAPTG', 'S4', ()), ('SPCRMPT', 'S4', ()),
             ('VRAMP', '<i8', ()), ('DEACTEL', 'S4', ()), ('IDACT', '<i8', ()), ('IDBUG', '<i8', ()),
             ('ENTHALP', 'S4', ()), ('CLUMP', '<i8', ()), ('ENLIN', '<i8', ()), ('MAPTOL', 'S4', ()),
             ('VMPTOL', '<f8', ()), ('INLAMC', 'S4', ()), ('VCOORD', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WAVELENGTH(object):
    name = 'WAVELENGTH'
    path = '/NASTRAN/INPUT/MATERIAL/RADBND'
    dtype = [('LAMBDA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATERIAL/RADBND'
    dtype = [('NUMBER', '<i8', ()), ('PLANCK2', '<f8', ()), ('POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATERIAL/RADBND/WAVELENGTH']


@register_table
class RADC(object):
    name = 'RADC'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('RADMID', '<i8', ()), ('EMIS', '<f8', ()), ('ABSO', '<f8', ()), ('IRSPEC', '<f8', ()),
             ('UVSPEC', '<f8', ()), ('FORMAT', 'S8', ()), ('IRTRANSP', '<f8', ()), ('IRTRANSL', '<f8', ()),
             ('UVTRANSP', '<f8', ()), ('UVTRANSL', '<f8', ()), ('IRREFI', '<f8', ()), ('UVREFI', '<f8', ()),
             ('ESTAR', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RADCT(object):
    name = 'RADCT'
    path = '/NASTRAN/INPUT/MATERIAL'
    dtype = [('RADMID', '<i8', ()), ('EMISID', '<i8', ()), ('ABSOID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EMIS(object):
    name = 'EMIS'
    path = '/NASTRAN/INPUT/MATERIAL/RADM'
    dtype = [('EMISI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class MIDS(object):
    name = 'MIDS'
    path = '/NASTRAN/INPUT/MATERIAL/RADM'
    dtype = [('RADMID', '<i8', ()), ('ABSORP', '<f8', ()), ('EMIS_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATERIAL/RADM'
    dtype = [('NUMBER', '<i8', ()), ('MIDS_POS', '<i8', ()), ('MIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATERIAL/RADM/EMIS', '/NASTRAN/INPUT/MATERIAL/RADM/MIDS']


@register_table
class MIDS(object):
    name = 'MIDS'
    path = '/NASTRAN/INPUT/MATERIAL/RADMT'
    dtype = [('RADMID', '<i8', ()), ('TABSORP', '<i8', ()), ('TEMIS_POS', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TEMIS(object):
    name = 'TEMIS'
    path = '/NASTRAN/INPUT/MATERIAL/RADMT'
    dtype = [('TA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATERIAL/RADMT'
    dtype = [('NUMBER', '<i8', ()), ('MIDS_POS', '<i8', ()), ('MIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATERIAL/RADMT/MIDS', '/NASTRAN/INPUT/MATERIAL/RADMT/TEMIS']


@register_table
class CONM1(object):
    name = 'CONM1'
    path = '/NASTRAN/INPUT/MATRIX'
    dtype = [('EID', '<i8', ()), ('G', '<i8', ()), ('CID', '<i8', ()), ('M1', '<f8', ()), ('M2', '<f8', (2,)),
             ('M3', '<f8', (3,)), ('M4', '<f8', (4,)), ('M5', '<f8', (5,)), ('M6', '<f8', (6,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJ'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJ'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJ'
    dtype = [('GJ', '<i8', ()), ('CJ', '<i8', ()), ('POS_I', '<i8', ()), ('NUM_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATRIX/DMIJ'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('NUM_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATRIX/DMIJ/ICLIST', '/NASTRAN/INPUT/MATRIX/DMIJ/IRLIST',
                 '/NASTRAN/INPUT/MATRIX/DMIJ/JLIST']


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJI'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJI'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIJI'
    dtype = [('GJ', '<i8', ()), ('CJ', '<i8', ()), ('POS_I', '<i8', ()), ('NUM_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATRIX/DMIJI'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('NUM_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATRIX/DMIJI/ICLIST', '/NASTRAN/INPUT/MATRIX/DMIJI/IRLIST',
                 '/NASTRAN/INPUT/MATRIX/DMIJI/JLIST']


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIK'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIK'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/MATRIX/DMIK'
    dtype = [('GJ', '<i8', ()), ('CJ', '<i8', ()), ('POS_I', '<i8', ()), ('NUM_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MATRIX/DMIK'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('NUM_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MATRIX/DMIK/ICLIST', '/NASTRAN/INPUT/MATRIX/DMIK/IRLIST',
                 '/NASTRAN/INPUT/MATRIX/DMIK/JLIST']


@register_table
class MASTERS(object):
    name = 'MASTERS'
    path = '/NASTRAN/INPUT/MODULES/MDBCNCT'
    dtype = [('MODMI', '<i8', ()), ('IDMA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SLAVES(object):
    name = 'SLAVES'
    path = '/NASTRAN/INPUT/MODULES/MDBCNCT'
    dtype = [('MODSI', '<i8', ()), ('IDSL', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDBCNCT'
    dtype = [('NWD', '<i8', ()), ('ID', '<i8', ()), ('BCGPID', '<i8', ()), ('BCPPID', '<i8', ()), ('MODS', '<i8', ()),
             ('IDSLAVE', '<i8', ()), ('MODM', '<i8', ()), ('IDMASTER', '<i8', ()), ('NSLA', '<i8', ()),
             ('NMAS', '<i8', ()), ('POS_SLA', '<i8', ()), ('POS_MAS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDBCNCT/MASTERS', '/NASTRAN/INPUT/MODULES/MDBCNCT/SLAVES']


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/MODULES/MDBCTB1'
    dtype = [('MODI', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDBCTB1'
    dtype = [('BCID', '<i8', ()), ('NID', '<i8', ()), ('IDLIST_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDBCTB1/IDLIST']


@register_table
class GIDA(object):
    name = 'GIDA'
    path = '/NASTRAN/INPUT/MODULES/MDBNDRY'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDBNDRY'
    dtype = [('MIDA', '<i8', ()), ('MIDB', '<i8', ()), ('GIDA_POS', '<i8', ()), ('GIDA_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDBNDRY/GIDA']


@register_table
class TOP(object):
    name = 'TOP'
    path = '/NASTRAN/INPUT/MODULES/MDBOLT'
    dtype = [('GTID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class BOTTOM(object):
    name = 'BOTTOM'
    path = '/NASTRAN/INPUT/MODULES/MDBOLT'
    dtype = [('GBID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDBOLT'
    dtype = [('ID', '<i8', ()), ('GRIDC', '<i8', ()), ('MODC', '<i8', ()), ('MODT', '<i8', ()), ('MODB', '<i8', ()),
             ('TOP_POS', '<i8', ()), ('TOP_LEN', '<i8', ()), ('BOTTOM_POS', '<i8', ()), ('BOTTOM_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDBOLT/TOP', '/NASTRAN/INPUT/MODULES/MDBOLT/BOTTOM']


@register_table
class MDBULK(object):
    name = 'MDBULK'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('MODID', '<i8', ()), ('METHOD', 'S8', ()), ('TOL', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/MODULES/MDCONCT'
    dtype = [('MID', '<i8', ()), ('GID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDCONCT'
    dtype = [('BID', '<i8', ()), ('TYPE', 'S8', ()), ('TOL', '<f8', ()), ('MODID', '<i8', ()), ('GRID', '<i8', ()),
             ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('CID', '<i8', ()), ('G_POS', '<i8', ()),
             ('G_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDCONCT/G']


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/MODULES/MDDMIG'
    dtype = [('MODJ', '<i8', ()), ('GJ', '<i8', ()), ('CJ', '<i8', ()), ('POS_I', '<i8', ()), ('LEN_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/MODULES/MDDMIG'
    dtype = [('MODI', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/MODULES/MDDMIG'
    dtype = [('MODI', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDDMIG'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('LEN_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDDMIG/JLIST', '/NASTRAN/INPUT/MODULES/MDDMIG/IRLIST',
                 '/NASTRAN/INPUT/MODULES/MDDMIG/ICLIST']


@register_table
class GID(object):
    name = 'GID'
    path = '/NASTRAN/INPUT/MODULES/MDEXCLD'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDEXCLD'
    dtype = [('MIDA', '<i8', ()), ('MIDB', '<i8', ()), ('GID_POS', '<i8', ()), ('GID_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDEXCLD/GID']


@register_table
class GCA(object):
    name = 'GCA'
    path = '/NASTRAN/INPUT/MODULES/MDMPC'
    dtype = [('MODI', '<i8', ()), ('GI', '<i8', ()), ('CI', '<i8', ()), ('AI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDMPC'
    dtype = [('SID', '<i8', ()), ('MOD1', '<i8', ()), ('G1', '<i8', ()), ('C1', '<i8', ()), ('A1', '<f8', ()),
             ('GCA_POS', '<i8', ()), ('GCA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDMPC/GCA']


@register_table
class MDRJNT(object):
    name = 'MDRJNT'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('EID', '<i8', ()), ('MODA', '<i8', ()), ('GA', '<i8', ()), ('MODB', '<i8', ()), ('GB', '<i8', ()),
             ('CB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MDFAST(object):
    name = 'MDFAST'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GS', '<i8', ()), ('FORMAT', 'S8', ()), ('GA', '<i8', ()),
             ('GB', '<i8', ()), ('MONE', '<i8', ()), ('IDA', '<i8', ()), ('IDB', '<i8', ()), ('XS', '<f8', ()),
             ('YS', '<f8', ()), ('ZS', '<f8', ()), ('MODA', '<i8', ()), ('MODB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MDLABEL(object):
    name = 'MDLABEL'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('MODID', '<i8', ()), ('LABEL', 'S56', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MODS(object):
    name = 'MODS'
    path = '/NASTRAN/INPUT/MODULES/MDRBE2'
    dtype = [('MODI', '<i8', ()), ('GMI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDRBE2'
    dtype = [('EID', '<i8', ()), ('MODN', '<i8', ()), ('GN', '<i8', ()), ('CM', '<i8', ()), ('ALPHA', '<f8', ()),
             ('MODS_POS', '<i8', ()), ('MODS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDRBE2/MODS']


@register_table
class WTCG(object):
    name = 'WTCG'
    path = '/NASTRAN/INPUT/MODULES/MDRBE3'
    dtype = [('WT', '<f8', ()), ('C', '<i8', ()), ('G_POS', '<i8', ()), ('G_LEN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class G(object):
    name = 'G'
    path = '/NASTRAN/INPUT/MODULES/MDRBE3'
    dtype = [('MODIJ', '<i8', ()), ('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class GM(object):
    name = 'GM'
    path = '/NASTRAN/INPUT/MODULES/MDRBE3'
    dtype = [('MODM', '<i8', ()), ('GM', '<i8', ()), ('CM', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/MODULES/MDRBE3'
    dtype = [('EID', '<i8', ()), ('REFMOD', '<i8', ()), ('REFG', '<i8', ()), ('REFC', '<i8', ()),
             ('WTCG_POS', '<i8', ()), ('WTCG_LEN', '<i8', ()), ('GM_POS', '<i8', ()), ('GM_LEN', '<i8', ()),
             ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/MODULES/MDRBE3/WTCG', '/NASTRAN/INPUT/MODULES/MDRBE3/G',
                 '/NASTRAN/INPUT/MODULES/MDRBE3/GM']


@register_table
class MDRROD(object):
    name = 'MDRROD'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('EID', '<i8', ()), ('MODA', '<i8', ()), ('GA', '<i8', ()), ('MODB', '<i8', ()), ('GB', '<i8', ()),
             ('CMA', '<i8', ()), ('CMB', '<i8', ()), ('ALPHA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MDSEAM(object):
    name = 'MDSEAM'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('SMLN', 'S8', ()), ('CTYPE', 'S8', ()), ('IDAS', '<i8', ()),
             ('IDBS', '<i8', ()), ('IDAE', '<i8', ()), ('IDBE', '<i8', ()), ('FLAG', '<i8', ()), ('GS', '<i8', ()),
             ('GE', '<i8', ()), ('XS', '<f8', ()), ('YS', '<f8', ()), ('ZS', '<f8', ()), ('XE', '<f8', ()),
             ('YE', '<f8', ()), ('ZE', '<f8', ()), ('MODA', '<i8', ()), ('MODB', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MDWELD(object):
    name = 'MDWELD'
    path = '/NASTRAN/INPUT/MODULES'
    dtype = [('EID', '<i8', ()), ('PID', '<i8', ()), ('GS', '<i8', ()), ('GA', '<i8', ()), ('GB', '<i8', ()),
             ('FORMAT', 'S8', ()), ('SPTYP', '<i8', ()), ('MCID', '<i8', ()), ('MODA', '<i8', ()), ('MODB', '<i8', ()),
             ('MODS', '<i8', ()), ('PIDA', '<i8', ()), ('PIDB', '<i8', ()), ('XS', '<f8', ()), ('YS', '<f8', ()),
             ('ZS', '<f8', ()), ('SHIDA', '<i8', ()), ('SHIDB', '<i8', ()), ('GIDA', '<i8', (8,)),
             ('GIDB', '<i8', (8,)), ('TAVG', '<f8', ()), ('TMIN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEGRID(object):
    name = 'AEGRID'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('GID', '<i8', ()), ('CP', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()),
             ('CD', '<i8', ()), ('SP', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EPOINT(object):
    name = 'EPOINT'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID(object):
    name = 'GRID'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('CP', '<i8', ()), ('X', '<f8', (3,)), ('CD', '<i8', ()), ('PS', '<i8', ()),
             ('SEID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDB(object):
    name = 'GRIDB'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('CD', '<i8', ()), ('PA', '<i8', ()), ('IDF', '<i8', ()), ('PHI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDF(object):
    name = 'GRIDF'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('R', '<f8', ()), ('Z', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('IDF', '<i8', ()), ('R', '<f8', ()), ('Z', '<f8', ()), ('W', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class POINTAX(object):
    name = 'POINTAX'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('RID', '<i8', ()), ('PHI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRESPT(object):
    name = 'PRESPT'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('IDF', '<i8', ()), ('IDP', '<i8', ()), ('PHI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RINGFL(object):
    name = 'RINGFL'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('IDF', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SEQGP(object):
    name = 'SEQGP'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('SEQID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPOINT(object):
    name = 'SPOINT'
    path = '/NASTRAN/INPUT/NODE'
    dtype = [('ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACMODL(object):
    name = 'ACMODL'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('INTER', 'S8', ()), ('INFOR', 'S8', ()), ('FSET', '<i8', ()), ('SSET', '<i8', ()), ('NORMAL', '<f8', ()),
             ('METHOD', 'S8', ()), ('SKNEPS', '<f8', ()), ('DSKNEPS', '<f8', ()), ('INTOL', '<f8', ()),
             ('ALLSSET', 'S8', ()), ('SRCHUNIT', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACPEMCP(object):
    name = 'ACPEMCP'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('TID', '<i8', ()), ('SGLUED', '<i8', ()), ('SSLIDE', '<i8', ()), ('SOPEN', '<i8', ()),
             ('SIMPER', '<i8', ()), ('OOC', '<i8', ()), ('SCUX', '<i8', ()), ('SCUY', '<i8', ()), ('SCUZ', '<i8', ()),
             ('SCRX', '<i8', ()), ('SCRY', '<i8', ()), ('SCRZ', '<i8', ()), ('SCFP', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LISTIDS(object):
    name = 'LISTIDS'
    path = '/NASTRAN/INPUT/PARAMETER/AECOMP'
    dtype = [('LISTID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/AECOMP'
    dtype = [('NAME', 'S8', ()), ('LISTTYPE', 'S8', ()), ('LISTIDS_POS', '<i8', ()), ('LISTIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/AECOMP/LISTIDS']


@register_table
class LABELS(object):
    name = 'LABELS'
    path = '/NASTRAN/INPUT/PARAMETER/AECOMPL'
    dtype = [('LABEL', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/AECOMPL'
    dtype = [('NAME', 'S8', ()), ('LABELS_POS', '<i8', ()), ('LABELS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/AECOMPL/LABELS']


@register_table
class DIS(object):
    name = 'DIS'
    path = '/NASTRAN/INPUT/PARAMETER/AEFACT'
    dtype = [('D', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/AEFACT'
    dtype = [('SID', '<i8', ()), ('DIS_POS', '<i8', ()), ('DIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/AEFACT/DIS']


@register_table
class AEPARM(object):
    name = 'AEPARM'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('UNIT', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AERO(object):
    name = 'AERO'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ACSID', '<i8', ()), ('VELOCITY', '<f8', ()), ('REFC', '<f8', ()), ('RHOREF', '<f8', ()),
             ('SYMXZ', '<i8', ()), ('SYMXY', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AEROS(object):
    name = 'AEROS'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ACSID', '<i8', ()), ('RCSID', '<i8', ()), ('REFC', '<f8', ()), ('REFB', '<f8', ()), ('REFS', '<f8', ()),
             ('SYMXZ', '<i8', ()), ('SYMXY', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AESCALE(object):
    name = 'AESCALE'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('X1REF', '<f8', ()), ('X2REF', '<f8', ()), ('X3REF', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AESTAT(object):
    name = 'AESTAT'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIC(object):
    name = 'AXIC'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('H', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HARMNRS(object):
    name = 'HARMNRS'
    path = '/NASTRAN/INPUT/PARAMETER/AXIF'
    dtype = [('N', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/AXIF'
    dtype = [('FCID', '<i8', ()), ('GRAV', '<f8', ()), ('DRHO', '<f8', ()), ('DBULK', '<f8', ()), ('NOSYM', '<i8', ()),
             ('HARMNRS_POS', '<i8', ()), ('HARMNRS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/AXIF/HARMNRS']


@register_table
class AXSLOT(object):
    name = 'AXSLOT'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('RHOD', '<f8', ()), ('BD', '<f8', ()), ('N', '<i8', ()), ('WD', '<f8', ()), ('MD', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TYPEFLAGS(object):
    name = 'TYPEFLAGS'
    path = '/NASTRAN/INPUT/PARAMETER/BLDOUT'
    dtype = [('ISLV', '<i8', ()), ('IMAST', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/BLDOUT'
    dtype = [('RID', '<i8', ()), ('KIND', '<i8', ()), ('TOUT', '<f8', ()), ('FILT', '<f8', ()), ('IFILE', '<i8', ()),
             ('MAP', '<i8', ()), ('TOL', '<f8', ()), ('MOM', '<i8', ()), ('MASS', '<f8', ()), ('ROFFSET', '<f8', ()),
             ('ZOFFSET', '<f8', ()), ('GID', '<i8', ()), ('ISET', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()),
             ('X3', '<f8', ()), ('TYPEFLAGS_POS', '<i8', ()), ('TYPEFLAGS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/BLDOUT/TYPEFLAGS']


@register_table
class MEMIDS(object):
    name = 'MEMIDS'
    path = '/NASTRAN/INPUT/PARAMETER/CASECC'
    dtype = [('SETMEM', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PARAMS(object):
    name = 'PARAMS'
    path = '/NASTRAN/INPUT/PARAMETER/CASECC'
    dtype = [('SID', '<i8', ()), ('NAME', 'S8', ()), ('TYPE', '<i8', ()), ('INT', '<i8', ()), ('CHAR', 'S8', ()),
             ('RE', '<f8', ()), ('IM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SETS(object):
    name = 'SETS'
    path = '/NASTRAN/INPUT/PARAMETER/CASECC'
    dtype = [('SID', '<i8', ()), ('SETID', '<i8', ()), ('SETLEN', '<i8', ()), ('MEMIDS_LEN', '<i8', ()),
             ('MEMIDS_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SUBCASE(object):
    name = 'SUBCASE'
    path = '/NASTRAN/INPUT/PARAMETER/CASECC'
    dtype = [('SID', '<i8', ()), ('MPCSET', '<i8', ()), ('SPCSET', '<i8', ()), ('ESLSET', '<i8', ()),
             ('REESET', '<i8', ()), ('ELDSET', '<i8', ()), ('THLDSET', '<i8', ()), ('THMATSET', '<i8', ()),
             ('TIC', '<i8', ()), ('NONPTSET', '<i8', ()), ('NONMEDIA', '<i8', ()), ('NONFMT', '<i8', ()),
             ('DYMLDSET', '<i8', ()), ('FEQRESET', '<i8', ()), ('TFSET', '<i8', ()), ('SYMFLG', '<i8', ()),
             ('LDSPTSET', '<i8', ()), ('LDSMEDIA', '<i8', ()), ('LDSFMT', '<i8', ()), ('DPLPTSET', '<i8', ()),
             ('DPLMEDIA', '<i8', ()), ('DPLFMT', '<i8', ()), ('STSPTSET', '<i8', ()), ('STSMEDIA', '<i8', ()),
             ('STSFMT', '<i8', ()), ('FCEPTSET', '<i8', ()), ('FCEMEDIA', '<i8', ()), ('FCEFMT', '<i8', ()),
             ('ACCPTSET', '<i8', ()), ('ACCMEDIA', '<i8', ()), ('ACCFMT', '<i8', ()), ('VELPTSET', '<i8', ()),
             ('VELMEDIA', '<i8', ()), ('VELFMT', '<i8', ()), ('FOCPTSET', '<i8', ()), ('FOCMEDIA', '<i8', ()),
             ('FOCFMT', '<i8', ()), ('TSTEPTRN', '<i8', ()), ('TITLE', 'S128', ()), ('SUBTITLE', 'S128', ()),
             ('LABEL', 'S128', ()), ('STPLTFLG', '<i8', ()), ('AXSYMSET', '<i8', ()), ('NOHARMON', '<i8', ()),
             ('TSTRU', '<i8', ()), ('SETKP', '<i8', ()), ('FLAGKP', '<i8', ()), ('K2PP', 'S8', ()),
             ('SETMP', '<i8', ()), ('FLAGMP', '<i8', ()), ('M2PP', 'S8', ()), ('SETBP', '<i8', ()),
             ('FLAGBP', '<i8', ()), ('B2PP', 'S8', ()), ('OUTRESPV', '<i8', ()), ('SEDR', '<i8', ()),
             ('FLDBNDY', '<i8', ()), ('CEESET', '<i8', ()), ('DAMPTBL', '<i8', ()), ('DYNRED', '<i8', ()),
             ('SSDSET', '<i8', ()), ('SSDMEDIA', '<i8', ()), ('SSDFMT', '<i8', ()), ('SSVSET', '<i8', ()),
             ('SSVMEDIA', '<i8', ()), ('SSVFMT', '<i8', ()), ('SSASET', '<i8', ()), ('SSAMEDIA', '<i8', ()),
             ('SSAFMT', '<i8', ()), ('NONLINLD', '<i8', ()), ('PARTIT', '<i8', ()), ('CYCLIC', '<i8', ()),
             ('RANDOM', '<i8', ()), ('NONPARAM', '<i8', ()), ('FLUTTER', '<i8', ()), ('LCC', '<i8', ()),
             ('GPFSET', '<i8', ()), ('GPFMEDIA', '<i8', ()), ('GPFFMT', '<i8', ()), ('ESESET', '<i8', ()),
             ('ESEMEDIA', '<i8', ()), ('ESEFMT', '<i8', ()), ('ARFPTSET', '<i8', ()), ('ARFMEDIA', '<i8', ()),
             ('ARFFMT', '<i8', ()), ('SEID', '<i8', ()), ('LCN', '<i8', ()), ('GUST', '<i8', ()),
             ('SEFINAL', '<i8', ()), ('SEMG', '<i8', ()), ('SEKR', '<i8', ()), ('SELG', '<i8', ()), ('SELR', '<i8', ()),
             ('SEEX', '<i8', ()), ('SETKG', '<i8', ()), ('FLAGKG', '<i8', ()), ('K2GG', 'S8', ()), ('SETMG', '<i8', ()),
             ('FLAGMG', '<i8', ()), ('M2GG', 'S8', ()), ('SETBG', '<i8', ()), ('FLAGBG', '<i8', ()), ('B2GG', 'S8', ()),
             ('SVSET', '<i8', ()), ('SVMEDIA', '<i8', ()), ('SVFMT', '<i8', ()), ('FLUPTSET', '<i8', ()),
             ('FLUMEDIA', '<i8', ()), ('FLUFMT', '<i8', ()), ('HOUT', '<i8', (3,)), ('NOUT', '<i8', (3,)),
             ('SETPG', '<i8', ()), ('FLAGPG', '<i8', ()), ('P2G', 'S8', ()), ('LOADSET', '<i8', ()),
             ('SEMR', '<i8', ()), ('VONMISES', '<i8', ()), ('SECMDFLG', '<i8', ()), ('GPSPTSET', '<i8', ()),
             ('GPSMEDIA', '<i8', ()), ('GPSFMT', '<i8', ()), ('STFSET', '<i8', ()), ('STFMEDIA', '<i8', ()),
             ('STFFMT', '<i8', ()), ('CLOAD', '<i8', ()), ('SET2ID', '<i8', ()), ('DSAPRT', '<i8', ()),
             ('DSASTORE', '<i8', ()), ('DSAOUTPT', '<i8', ()), ('STNSET', '<i8', ()), ('STNMEDIA', '<i8', ()),
             ('STNFMT', '<i8', ()), ('APRESS', '<i8', ()), ('TRIM', '<i8', ()), ('MODLIST', '<i8', ()),
             ('REESETF', '<i8', ()), ('ESDPTSET', '<i8', ()), ('ESDMEDIA', '<i8', ()), ('ESDFMT', '<i8', ()),
             ('GSDPTSET', '<i8', ()), ('GSDMEDIA', '<i8', ()), ('GSDFMT', '<i8', ()), ('SEDV', '<i8', ()),
             ('SERE', '<i8', ()), ('SERS', '<i8', ()), ('CNTSET', '<i8', ()), ('CNTMEDIA', '<i8', ()),
             ('CNTFMT', '<i8', ()), ('DIVERG', '<i8', ()), ('OUTRCV', '<i8', ()), ('STATSUBP', '<i8', ()),
             ('DFT1', '<f8', ()), ('DFT2', '<f8', ()), ('ADAPT', '<i8', ()), ('DESOBJ', '<i8', ()),
             ('DESSUB', '<i8', ()), ('SUBSPAN', '<i8', ()), ('DESGLB', '<i8', ()), ('ANALYSIS', 'S4', ()),
             ('GPQSTRS', '<i8', ()), ('GPQFORC', '<i8', ()), ('GPQSTRN', '<i8', ()), ('SUPORT1', '<i8', ()),
             ('STATSUBB', '<i8', ()), ('BCID', '<i8', ()), ('AUXMODEL', '<i8', ()), ('ADACT', '<i8', ()),
             ('DATSET', '<i8', ()), ('DATMEDIA', '<i8', ()), ('DATFMT', '<i8', ()), ('VUGSET', '<i8', ()),
             ('VUGMEDIA', '<i8', ()), ('VUGFMT', '<i8', ()), ('MPCFSET', '<i8', ()), ('MPCMEDIA', '<i8', ()),
             ('MPCFFMT', '<i8', ()), ('REUESET', '<i8', ()), ('DAMPTBLF', '<i8', ()), ('ITERMETH', '<i8', ()),
             ('NLSSET', '<i8', ()), ('NLSMEDIA', '<i8', ()), ('NLSFMT', '<i8', ()), ('MODTRKID', '<i8', ()),
             ('DSAFORM', '<i8', ()), ('DSAEXPO', '<i8', ()), ('DSABEGIN', '<i8', ()), ('DSAINTVL', '<i8', ()),
             ('DSAFINAL', '<i8', ()), ('DSASETID', '<i8', ()), ('SORTFLG', '<i8', ()), ('RANDBIT', '<i8', ()),
             ('AECONFIG', 'S8', ()), ('AESYMXY', '<i8', ()), ('AESYMXZ', '<i8', ()), ('OCIDREQ', '<i8', ()),
             ('GPEPTSET', '<i8', ()), ('GPEMEDIA', '<i8', ()), ('GPEFMT', '<i8', ()), ('TEMPMAT', '<i8', ()),
             ('AECSSSET', '<i8', ()), ('EKEPTSET', '<i8', ()), ('EKEMEDIA', '<i8', ()), ('EKEFMT', '<i8', ()),
             ('EKETHRSH', '<f8', ()), ('EDEPTSET', '<i8', ()), ('EDEMEDIA', '<i8', ()), ('EDEFMT', '<i8', ()),
             ('EDETHRSH', '<f8', ()), ('DFT3', '<f8', ()), ('DFR1', '<f8', ()), ('DFR2', '<f8', ()),
             ('DFR3', '<f8', ()), ('SETK4G', '<i8', ()), ('FLAGK4G', '<i8', ()), ('K42GG', 'S8', ()),
             ('A2GG', 'S8', ()), ('NK42GG', '<i8', ()), ('NA2GG', '<i8', ()), ('EFFMASET', '<i8', ()),
             ('EFFMAGID', '<i8', ()), ('EFFMATHR', '<f8', ()), ('EQUILMED', '<i8', ()), ('EQUILGRD', '<i8', ()),
             ('RCRSET', '<i8', ()), ('RCRFMT', '<i8', ()), ('AEUXREF', '<i8', ()), ('GCHK', '<i8', ()),
             ('GCHKOUT', '<i8', ()), ('GCHKSET', '<i8', ()), ('GCHKGID', '<i8', ()), ('GCHKTHR', '<f8', ()),
             ('GCHKRTHR', '<f8', ()), ('GCHKDREC', '<i8', ()), ('ASPCMED', '<i8', ()), ('ASPCEPS', '<f8', ()),
             ('ASPCPRT', '<f8', ()), ('ASPCPCH', '<i8', ()), ('NK2PP', '<i8', ()), ('NM2PP', '<i8', ()),
             ('NB2PP', '<i8', ()), ('NK2GG', '<i8', ()), ('NM2GG', '<i8', ()), ('NB2GG', '<i8', ()),
             ('NP2G', '<i8', ()), ('GEODSET', '<i8', ()), ('GEODMXMN', '<i8', ()), ('GEODOCID', '<i8', ()),
             ('GEODNUMB', '<i8', ()), ('GEOLSET', '<i8', ()), ('GEOLMXMN', '<i8', ()), ('GEOLOCID', '<i8', ()),
             ('GEOLNUMB', '<i8', ()), ('GEOSSET', '<i8', ()), ('GEOSMXMN', '<i8', ()), ('GEOSOCID', '<i8', ()),
             ('GEOSNUMB', '<i8', ()), ('GEOMSET', '<i8', ()), ('GEOMMXMN', '<i8', ()), ('GEOMOCID', '<i8', ()),
             ('GEOMNUMB', '<i8', ()), ('GEOASET', '<i8', ()), ('GEOAMXMN', '<i8', ()), ('GEOAOCID', '<i8', ()),
             ('GEOANUMB', '<i8', ()), ('GEOVSET', '<i8', ()), ('GEOVMXMN', '<i8', ()), ('GEOVOCID', '<i8', ()),
             ('GEOVNUMB', '<i8', ()), ('NTFL', '<i8', ()), ('GPKESET', '<i8', ()), ('GPKEMEDI', '<i8', ()),
             ('GPKEFMT', '<i8', ()), ('SEDAMP', '<i8', ()), ('WCHK', '<i8', ()), ('WCHKOUT', '<i8', ()),
             ('WCHKSET', '<i8', ()), ('WCHKGID', '<i8', ()), ('WCHKCGI', '<i8', ()), ('WCHKWM', '<i8', ()),
             ('EXSEOUT', '<i8', ()), ('EXSEMED', '<i8', ()), ('EXSEUNIT', '<i8', ()), ('EXSERES1', '<i8', ()),
             ('EXSERES2', '<i8', ()), ('FK2PP', '<i8', ()), ('FM2PP', '<i8', ()), ('FB2PP', '<i8', ()),
             ('FK2GG', '<i8', ()), ('FM2GG', '<i8', ()), ('FB2GG', '<i8', ()), ('TICTYPE', 'S4', ()),
             ('FK42GG', '<i8', ()), ('FA2GG', '<i8', ()), ('SUBSTEP', '<i8', ()), ('STEPID', '<i8', ()),
             ('NSMID', '<i8', ()), ('ROUTDISP', '<i8', ()), ('ROUTVELO', '<i8', ()), ('ROUTACCE', '<i8', ()),
             ('ROUTLOAD', '<i8', ()), ('ROUTSPCF', '<i8', ()), ('ROUTSTRS', '<i8', ()), ('ROUTFORC', '<i8', ()),
             ('ROUTSTRN', '<i8', ()), ('ROUTMSCF', '<i8', ()), ('MDLSSET', '<i8', ()), ('MDLSMEDIA', '<i8', ()),
             ('MDLSFMT', '<i8', ()), ('MDLSESRT', 'S4', ()), ('MDLSTHRE', '<f8', ()), ('MDLSTFVL', '<i8', ()),
             ('MDLKSET', '<i8', ()), ('MDLKMEDIA', '<i8', ()), ('MDLKFMT', '<i8', ()), ('MDLKESRT', 'S4', ()),
             ('MDLKTHRE', '<f8', ()), ('MDLKTFVL', '<i8', ()), ('ACPOWSET', '<i8', ()), ('ACPOWMED', '<i8', ()),
             ('ACPOWFMT', '<i8', ()), ('ACPOWCSV', '<i8', ()), ('NLOUT', '<i8', ()), ('SEEFMNO', '<i8', ()),
             ('SEEFMHV', 'S4', ()), ('SEEFMDLF', '<i8', ()), ('SEEFMBND', '<i8', ()), ('SEEFMSDP', '<i8', ()),
             ('CONNECTOR', '<i8', ()), ('ESETHRSH', '<f8', ()), ('POSTUNIT', '<i8', ()), ('POSTOPT1', '<i8', ()),
             ('POSTOPT2', '<i8', ()), ('TICDIFF', 'S4', ()), ('DSAESEID', '<i8', ()), ('MXMNGSET', '<i8', ()),
             ('MXMNGMDA', '<i8', ()), ('MXMNGFMT', '<i8', ()), ('MXMNESET', '<i8', ()), ('MXMNEMDA', '<i8', ()),
             ('MXMNEFMT', '<i8', ()), ('MCFRSET', '<i8', ()), ('MCFRSOLN', '<i8', ()), ('MCFRFILT', '<f8', ()),
             ('MCFROPT', '<i8', ()), ('ELSUMID', '<i8', ()), ('ELSUMOPT', '<i8', ()), ('ELSUMDUM', '<i8', ()),
             ('RGYRO', '<i8', ()), ('CMSESET', '<i8', ()), ('CMSEMDIA', '<i8', ()), ('CMSEOPTS', '<i8', ()),
             ('CMSETHRE', '<f8', ()), ('CMSETOPN', '<i8', ()), ('GPRSORT', '<i8', ()), ('MASSSET', '<i8', ()),
             ('AESOLN', 'S8', ()), ('POSTO2NM', 'S8', ()), ('RANDVAR', '<i8', ()), ('RSVCRQTS', '<i8', ()),
             ('RSVCOPTS', '<i8', ()), ('RSVCSTBS', '<i8', ()), ('RSVCRQTC', '<i8', ()), ('RSVCOPTC', '<i8', ()),
             ('RSVCSTBC', '<i8', ()), ('DESVAR', '<i8', ()), ('BCONTACTI', '<i8', ()), ('BCONTACTC', 'S4', ()),
             ('MODSELS1', '<i8', ()), ('MODSELS2', '<i8', ()), ('MODSELS3', '<i8', ()), ('MODSELS4', '<f8', ()),
             ('MODSELS5', '<f8', ()), ('MODSELS6', '<f8', ()), ('MODSELS7', '<f8', ()), ('MODSELS8', '<f8', ()),
             ('MODSELS9', '<f8', ()), ('MODSELF1', '<i8', ()), ('MODSELF2', '<i8', ()), ('MODSELF3', '<i8', ()),
             ('MODSELF4', '<f8', ()), ('MODSELF5', '<f8', ()), ('MODSELF6', '<f8', ()), ('MODSELF7', '<f8', ()),
             ('MODSELF8', '<f8', ()), ('MODSELF9', '<f8', ()), ('FTNURN', '<i8', ()), ('SUFNAM1', '<i8', ()),
             ('SUFNAM2', '<i8', ()), ('ENVELOP1', 'S4', ()), ('ENVELOP2', 'S4', ()), ('GPFLXSET', '<i8', ()),
             ('GPFLXMED', '<i8', ()), ('CAMPBELL', '<i8', ()), ('SPLINOUT', '<i8', (2,)), ('MONITOR', '<i8', ()),
             ('FBODYLD', '<i8', (2,)), ('STOCHAST', '<i8', ()), ('EXPTLDID', '<i8', ()), ('EXPTLDNM', 'S8', ()),
             ('EXPTLDSI', '<i8', ()), ('AERCONFIG', 'S8', ()), ('NLHARM', '<i8', ()), ('PFMSSID', '<i8', ()),
             ('PFMSMED', '<i8', ()), ('PFMSFMT', '<i8', ()), ('PFMSSSET', '<i8', ()), ('PFMSFLTR', '<f8', ()),
             ('PFMSOPTS', '<i8', ()), ('PFMSSTMP', '<i8', ()), ('PFMFSID', '<i8', ()), ('PFMFMED', '<i8', ()),
             ('PFMFFMT', '<i8', ()), ('PFMFSSET', '<i8', ()), ('PFMFPSET', '<i8', ()), ('PFMFFLTR', '<f8', ()),
             ('PFMFOPTS', '<i8', ()), ('PFMFFLMP', '<i8', ()), ('PFMFSTMP', '<i8', ()), ('PFPSID', '<i8', ()),
             ('PFPMED', '<i8', ()), ('PFPFMT', '<i8', ()), ('PFPSSET', '<i8', ()), ('PFPFLTR', '<f8', ()),
             ('PFPOPTS', '<i8', ()), ('PFPPSET', '<i8', ()), ('PFGSID', '<i8', ()), ('PFGMED', '<i8', ()),
             ('PFGFMT', '<i8', ()), ('PFGSSET', '<i8', ()), ('PFGGSET', '<i8', ()), ('NLICCASE', '<i8', ()),
             ('NLICSTEP', '<i8', ()), ('NLICLFAC', '<f8', ()), ('ACFPMSET', '<i8', ()), ('ACFPMMED', '<i8', ()),
             ('ACFPMFMT', '<i8', ()), ('FRFFLAG', '<i8', ()), ('FRFCMPID', '<i8', ()), ('FRFCONST', '<i8', ()),
             ('FRFUNTNO', '<i8', ()), ('FRFCMPNM', '<i8', (2,)), ('DESMOD', 'S8', ()), ('RSDAMPST', '<i8', ()),
             ('RSDAMPFL', '<i8', ()), ('VCCT', '<i8', ()), ('FP2G', '<i8', ()), ('FRQVAR', '<i8', ()),
             ('HADAPT', '<i8', ()), ('BCHANGE', '<i8', ()), ('BCMOVE', '<i8', ()), ('BSQUEAL', '<i8', ()),
             ('UNGLUE', '<i8', ()), ('HSUBCASE', '<i8', ()), ('HSTEP', '<i8', ()), ('HTIME', '<f8', ()),
             ('ERPSID', '<i8', ()), ('ERPMED', '<i8', ()), ('ERPFMT', '<i8', ()), ('ERPSSET', '<i8', ()),
             ('ERPFLTR', '<f8', ()), ('ERPOPTS', '<i8', ()), ('ERPCSV', '<i8', ()), ('TESTTHRR', '<f8', ()),
             ('TESTTHRI', '<i8', ()), ('ASMOUTFL', '<i8', ()), ('ASMOUTNM', '<i8', (2,)), ('ICFUNTNO', '<i8', ()),
             ('TFOSET', '<i8', ()), ('TFUNIT', '<i8', ()), ('TFLSET', '<i8', ()), ('TFOPTS', '<i8', ()),
             ('NLOOPH', '<i8', ()), ('NLSTEP', '<i8', ()), ('NMODES', '<i8', ()), ('RCPARM', '<i8', ()),
             ('NLOPCTRL', '<i8', ()), ('NLOPDBG', '<i8', ()), ('NLOPPOST', '<i8', ()), ('NLOPMPCH', '<i8', ()),
             ('NLICSSTP', '<i8', ()), ('DEACTEL', '<i8', ()), ('ACTIVAT', '<i8', ()), ('INTENSET', '<i8', ()),
             ('INTENMED', '<i8', ()), ('INTENFMT', '<i8', ()), ('NLICTOLR', '<f8', ()), ('PFPSSID', '<i8', ()),
             ('PFPSMED', '<i8', ()), ('PFPSFMT', '<i8', ()), ('PFPSSSET', '<i8', ()), ('PFPSFLTR', '<f8', ()),
             ('PFPSOPTS', '<i8', ()), ('PFPSPSET', '<i8', ()), ('ACTISET', '<i8', ()), ('ELSOSET', '<i8', ()),
             ('ELSRSET', '<i8', ()), ('ELSSSET', '<i8', ()), ('ELSTHRS', '<f8', ()), ('ELSBITS', '<i8', ()),
             ('WTSOSET', '<i8', ()), ('WTSRSET', '<i8', ()), ('WTSSSET', '<i8', ()), ('WTSTHRS', '<f8', ()),
             ('WTSBITS', '<i8', ()), ('PACOSET', '<i8', ()), ('PACSSET', '<i8', ()), ('PACBITS', '<i8', ()),
             ('IRLOAD', '<i8', ()), ('ICFSET', '<i8', ()), ('ICFMED', '<i8', ()), ('ICFFMT', '<i8', ()),
             ('ICFGENST', '<i8', ()), ('ICFGENNM', '<i8', (2,)), ('ICFUSEST', '<i8', ()), ('ICFUSENM', '<i8', (2,)),
             ('HISTSET', '<i8', ()), ('HISTTYPE', '<i8', ()), ('HISTFMT', '<i8', ()), ('FATIGUE', '<i8', ()),
             ('FTGMED', '<i8', ()), ('FTGFMT', '<i8', ()), ('NLOPDELI', '<i8', ()), ('NLOPGRID', '<i8', ()),
             ('DASAVE', '<i8', ()), ('GVSET', '<i8', ()), ('GVMEDIA', '<i8', ()), ('GVFMT', '<i8', ()),
             ('ERMPF', '<i8', ()), ('NVELOSET', '<i8', ()), ('NVELFDEF', '<i8', ()), ('NVELFDTA', '<i8', ()),
             ('NVELTHRS', '<f8', ()), ('NVELBITS', '<i8', ()), ('VITSOSET', '<i8', ()), ('VITSRSET', '<i8', ()),
             ('VITSSSET', '<i8', ()), ('VITSTHRS', '<f8', ()), ('VITSBITS', '<i8', ()), ('THLDVER', '<i8', ()),
             ('THMATVER', '<i8', ()), ('NPEAK', '<i8', ()), ('NEAR', '<f8', ()), ('LFREQ', '<f8', ()),
             ('HFREQ', '<f8', ()), ('RTYPE', '<i8', ()), ('PSCALE', '<i8', ()), ('NSAMP', '<i8', ()),
             ('MONSET', '<i8', ()), ('SEED', '<i8', ()), ('OFFD', '<f8', ()), ('FSORT2', '<i8', ()),
             ('LSEM', '<i8', ()), ('SYM_LEN', '<i8', ()), ('SYM_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SYM(object):
    name = 'SYM'
    path = '/NASTRAN/INPUT/PARAMETER/CASECC'
    dtype = [('COEF', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CBARAO(object):
    name = 'CBARAO'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('EID', '<i8', ()), ('SCALE', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()),
             ('X4', '<f8', ()), ('X5', '<f8', ()), ('X6', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DAMPING(object):
    name = 'DAMPING'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('G', '<f8', ()), ('ALPHA1', '<f8', ()), ('ALPHA2', '<f8', ()), ('ALPHA3', '<f8', ()),
             ('HYBRID', '<i8', ()), ('GEFACT', '<f8', ()), ('ROTSEP', 'S4', ()), ('W3', '<f8', ()), ('W4', '<f8', ()),
             ('WH', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MACHS(object):
    name = 'MACHS'
    path = '/NASTRAN/INPUT/PARAMETER/DIVERG'
    dtype = [('MI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/DIVERG'
    dtype = [('SID', '<i8', ()), ('NROOT', '<i8', ()), ('MACHS_POS', '<i8', ()), ('MACHS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/DIVERG/MACHS']


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIAX'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('NI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIAX'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('NI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIAX'
    dtype = [('GJ', '<i8', ()), ('CJ', '<i8', ()), ('NJ', '<i8', ()), ('POS_I', '<i8', ()), ('NUM_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/DMIAX'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('NUM_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/DMIAX/ICLIST', '/NASTRAN/INPUT/PARAMETER/DMIAX/IRLIST',
                 '/NASTRAN/INPUT/PARAMETER/DMIAX/JLIST']


@register_table
class ICLIST(object):
    name = 'ICLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIG'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VCR', '<f8', ()), ('VCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IRLIST(object):
    name = 'IRLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIG'
    dtype = [('GI', '<i8', ()), ('CI', '<i8', ()), ('VR', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class JLIST(object):
    name = 'JLIST'
    path = '/NASTRAN/INPUT/PARAMETER/DMIG'
    dtype = [('GJ', '<i8', ()), ('CJ', '<i8', ()), ('POS_I', '<i8', ()), ('NUM_I', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/DMIG'
    dtype = [('NAME', 'S8', ()), ('MATFORM', '<i8', ()), ('MATTYPE', '<i8', ()), ('MATCOLS', '<i8', ()),
             ('POS_J', '<i8', ()), ('NUM_J', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/DMIG/ICLIST', '/NASTRAN/INPUT/PARAMETER/DMIG/IRLIST',
                 '/NASTRAN/INPUT/PARAMETER/DMIG/JLIST']


@register_table
class DOPTPRM2(object):
    name = 'DOPTPRM2'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('APRCOD', '<i8', ()), ('IPRINT', '<i8', ()), ('METHOD', '<i8', ()), ('DELP', '<f8', ()),
             ('DPMIN', '<f8', ()), ('PTOL', '<f8', ()), ('CONV1', '<f8', ()), ('CONV2', '<f8', ()), ('GMAX', '<f8', ()),
             ('DELB', '<f8', ()), ('DELBM', '<f8', ()), ('CONVDV', '<f8', ()), ('CONVPR', '<f8', ()), ('P1', '<i8', ()),
             ('P2', '<i8', ()), ('DELX', '<f8', ()), ('DXMIN', '<f8', ()), ('ISCAL', '<i8', ()), ('CT', '<f8', ()),
             ('CTMIN', '<f8', ()), ('IPRNT1', '<i8', ()), ('IPRNT2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FIS(object):
    name = 'FIS'
    path = '/NASTRAN/INPUT/PARAMETER/FLFACT'
    dtype = [('FI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/FLFACT'
    dtype = [('SID', '<i8', ()), ('FIS_POS', '<i8', ()), ('FIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/FLFACT/FIS']


@register_table
class FLSYM(object):
    name = 'FLSYM'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('M', '<i8', ()), ('S1', '<i8', ()), ('S2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLUTTER(object):
    name = 'FLUTTER'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('METHOD', 'S8', ()), ('DENS', '<i8', ()), ('MACH', '<i8', ()), ('RFREQ', '<i8', ()),
             ('IMETH', 'S8', ()), ('SFLG', '<i8', ()), ('NVALUE', '<i8', ()), ('OMAX', '<f8', ()), ('EPS', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FSICTRL(object):
    name = 'FSICTRL'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SERV_ID', 'S8', ()), ('TYPE', 'S8', ()), ('FREQ', '<i8', ()), ('ANALYSIS', 'S8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HADACRI(object):
    name = 'HADACRI'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('CRITID', '<i8', ()), ('TYPE', '<i8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()), ('F3', '<f8', ()),
             ('F4', '<f8', ()), ('F5', '<f8', ()), ('F6', '<f8', ()), ('F7', '<f8', ()), ('F8', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HADAPTG(object):
    name = 'HADAPTG'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('FREQ', '<i8', ()), ('CONT', 'S8', ()), ('REPEAT', '<i8', ()), ('CRITID', '<i8', ()),
             ('WHEREMET', 'S4', ()), ('WHEREID', '<i8', ()), ('MSHGENID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HADAPTL(object):
    name = 'HADAPTL'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('FREQ', '<i8', ()), ('CONT', 'S8', ()), ('REPEAT', '<i8', ()), ('CRITID', '<i8', ()),
             ('WHEREMET', 'S4', ()), ('WHEREID', '<i8', ()), ('SNAPMTHD', '<i8', ()), ('MAXLEVEL', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HYBDAMP(object):
    name = 'HYBDAMP'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('METHOD', '<i8', ()), ('SDAMP', '<i8', ()), ('KDAMP', 'S4', ()), ('PRTEIG', 'S4', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class STRAINS(object):
    name = 'STRAINS'
    path = '/NASTRAN/INPUT/PARAMETER/IPSTRAIN'
    dtype = [('STRAIN', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/IPSTRAIN'
    dtype = [('EID1', '<i8', ()), ('EID2', '<i8', ()), ('INT1', '<i8', ()), ('INTN', '<i8', ()), ('LAY1', '<i8', ()),
             ('LAYN', '<i8', ()), ('STRAINS_POS', '<i8', ()), ('STRAINS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/IPSTRAIN/STRAINS']


@register_table
class STRESSES(object):
    name = 'STRESSES'
    path = '/NASTRAN/INPUT/PARAMETER/ISTRESS'
    dtype = [('STRESS', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/ISTRESS'
    dtype = [('EID1', '<i8', ()), ('EID2', '<i8', ()), ('INT1', '<i8', ()), ('INTN', '<i8', ()), ('LAY1', '<i8', ()),
             ('LAYN', '<i8', ()), ('STRESSES_POS', '<i8', ()), ('STRESSES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/ISTRESS/STRESSES']


@register_table
class ITER(object):
    name = 'ITER'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('PRECOND', 'S8', ()), ('CONV', 'S4', ()), ('MSGFLG', 'S4', ()), ('ITSEPS', '<f8', ()),
             ('ITSMAX', '<i8', ()), ('IPAD', '<i8', ()), ('IEXT', '<i8', ()), ('PREFONLY', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MDLPRM(object):
    name = 'MDLPRM'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('NAME', 'S8', ()), ('VALUE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CAERIDS(object):
    name = 'CAERIDS'
    path = '/NASTRAN/INPUT/PARAMETER/MONCNCM'
    dtype = [('CAERID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/MONCNCM'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('MREF', '<f8', ()), ('CAERIDS_POS', '<i8', ()),
             ('CAERIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/MONCNCM/CAERIDS']


@register_table
class MONDSP1(object):
    name = 'MONDSP1'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('AXES', '<i8', ()), ('COMP', 'S8', ()), ('CP', '<i8', ()),
             ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('CD', '<i8', ()), ('INDDOF', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MONPNTS(object):
    name = 'MONPNTS'
    path = '/NASTRAN/INPUT/PARAMETER/MONGRP'
    dtype = [('NAME', 'S8', ()), ('CLASS', 'S8', ()), ('RVAL', '<f8', ()), ('IVAL', '<i8', ()), ('STRNG', 'S32', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/MONGRP'
    dtype = [('GRPNAM', 'S8', ()), ('GRPLABEL', 'S56', ()), ('MONPNTS_POS', '<i8', ()), ('MONPNTS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/MONGRP/MONPNTS']


@register_table
class MONPNT1(object):
    name = 'MONPNT1'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('AXES', '<i8', ()), ('COMP', 'S8', ()), ('CP', '<i8', ()),
             ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('CD', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MONPNT2(object):
    name = 'MONPNT2'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('TABLE', 'S8', ()), ('TYPE', 'S8', ()), ('ITEM', 'S8', ()),
             ('EID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MONPNT3(object):
    name = 'MONPNT3'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('AXES', '<i8', ()), ('GRIDSET', '<i8', ()),
             ('ELEMSET', '<i8', ()), ('CID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()),
             ('XFLAG', '<i8', ()), ('CD', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COMPONENT(object):
    name = 'COMPONENT'
    path = '/NASTRAN/INPUT/PARAMETER/MONSUM'
    dtype = [('AXES', '<i8', ()), ('TYPE_POS', '<i8', ()), ('TYPE_LEN', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TYPE(object):
    name = 'TYPE'
    path = '/NASTRAN/INPUT/PARAMETER/MONSUM'
    dtype = [('CLASS', 'S8', ()), ('NAME', 'S8', ()), ('AXES', '<i8', ()), ('COEF', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/MONSUM'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('CP', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()),
             ('Z', '<f8', ()), ('CD', '<i8', ()), ('CLASS', 'S8', ()), ('COMPONENT_POS', '<i8', ()),
             ('COMPONENT_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/MONSUM/COMPONENT', '/NASTRAN/INPUT/PARAMETER/MONSUM/TYPE']


@register_table
class TYPE(object):
    name = 'TYPE'
    path = '/NASTRAN/INPUT/PARAMETER/MONSUMT'
    dtype = [('CLASS', 'S8', ()), ('NAME1', 'S8', ()), ('NAME2', 'S8', ()), ('NAME3', 'S8', ()), ('NAME4', 'S8', ()),
             ('NAME5', 'S8', ()), ('NAME6', 'S8', ()), ('NAME7', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/MONSUMT'
    dtype = [('NAME', 'S8', ()), ('LABEL', 'S56', ()), ('AXES', '<i8', ()), ('CP', '<i8', ()), ('X', '<f8', ()),
             ('Y', '<f8', ()), ('Z', '<f8', ()), ('CD', '<i8', ()), ('TYPE_POS', '<i8', ()), ('TYPE_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/MONSUMT/TYPE']


@register_table
class NLPARM(object):
    name = 'NLPARM'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('NINC', '<i8', ()), ('DT', '<f8', ()), ('KMETHOD', '<i8', ()), ('KSTEP', '<i8', ()),
             ('MAXITER', '<i8', ()), ('CONV', '<i8', ()), ('INTOUT', '<i8', ()), ('EPSU', '<f8', ()),
             ('EPSP', '<f8', ()), ('EPSW', '<f8', ()), ('MAXDIV', '<i8', ()), ('MAXQN', '<i8', ()),
             ('MAXLS', '<i8', ()), ('FSTRESS', '<f8', ()), ('LSTOL', '<f8', ()), ('MAXBIS', '<i8', ()),
             ('MAXR', '<f8', ()), ('RTOLB', '<f8', ()), ('MINITER', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NLPCI(object):
    name = 'NLPCI'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S4', ()), ('MINALR', '<f8', ()), ('MAXALR', '<f8', ()), ('SCALE', '<f8', ()),
             ('DESITER', '<i8', ()), ('MXINC', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class NLSTEP(object):
    name = 'NLSTEP'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('ID', '<i8', ()), ('TOTTIME', '<f8', ()), ('CTRLDEF', 'S8', ()), ('GENERAL', 'S8', ()),
             ('NBMWDSG', '<i8', ()), ('MAXITER', '<i8', ()), ('MINITER', '<i8', ()), ('MAXBIS', '<i8', ()),
             ('CREEP', '<i8', ()), ('FIXED', 'S8', ()), ('NBMWDSF', '<i8', ()), ('NINC', '<i8', ()), ('NO', '<i8', ()),
             ('ADAPT', 'S8', ()), ('NBMWDSA', '<i8', ()), ('DTINTF', '<f8', ()), ('DTMINF', '<f8', ()),
             ('DTMAXF', '<f8', ()), ('NDESIR', '<i8', ()), ('SFACT', '<f8', ()), ('INTOUT', '<i8', ()),
             ('NSMAX', '<i8', ()), ('IDAMP', '<i8', ()), ('DAMP', '<f8', ()), ('CRITTID', '<i8', ()),
             ('IPHYS', '<i8', ()), ('LIMTAR', '<i8', ()), ('RSMALL', '<f8', ()), ('RBIG', '<f8', ()),
             ('ADJUST', '<i8', ()), ('MSTEP', '<i8', ()), ('RB', '<f8', ()), ('UTOL', '<f8', ()), ('MECH', 'S8', ()),
             ('NBMWDSM', '<i8', ()), ('CONV', '<i8', ()), ('EPSU', '<f8', ()), ('EPSP', '<f8', ()), ('EPSW', '<f8', ()),
             ('KMETHOD', 'S4', ()), ('KSTEP', '<i8', ()), ('MRCONV', '<i8', ()), ('MAXQN', '<i8', ()),
             ('MAXLS', '<i8', ()), ('LSTOL', '<f8', ()), ('FSTRESS', '<f8', ()), ('HEAT', 'S8', ()),
             ('NBMWDSH', '<i8', ()), ('CONVH', '<i8', ()), ('EPSUH', '<f8', ()), ('EPSPH', '<f8', ()),
             ('EPSWH', '<f8', ()), ('KMETHODH', 'S4', ()), ('KSTEPH', '<i8', ()), ('MAXQNH', '<i8', ()),
             ('MAXLSH', '<i8', ()), ('LSTOLH', '<f8', ()), ('COUP', 'S8', ()), ('NBMWDSC', '<i8', ()),
             ('HGENPLAS', '<f8', ()), ('HGENFRIC', '<f8', ()), ('ARCLN', 'S8', ()), ('NBMWDSAR', '<i8', ()),
             ('TYPE', 'S8', ()), ('DTINITFA', '<f8', ()), ('MINALR', '<f8', ()), ('MAXALR', '<f8', ()),
             ('SCALEA', '<f8', ()), ('NDESIRA', '<i8', ()), ('NSMAXA', '<i8', ()), ('RCHEAT', 'S8', ()),
             ('NBMWDSR', '<i8', ()), ('SOLVER', 'S8', ()), ('DRLXCA', '<f8', ()), ('ARLXCA', '<f8', ()),
             ('BALENG', '<f8', ()), ('DAMPC', '<f8', ()), ('GRVCON', '<f8', ()), ('CSGFAC', '<f8', ()),
             ('NRLOOP', '<i8', ()), ('OUTINT', '<f8', ()), ('DTIMEI', '<f8', ()), ('LCNT', 'S8', ()),
             ('NBMWDSL', '<i8', ()), ('NINCC', '<i8', ()), ('CONVC', '<i8', ()), ('EPSUC', '<f8', ()),
             ('EPSPC', '<f8', ()), ('EPSWC', '<f8', ()), ('MAXDIVC', '<i8', ()), ('MAXBISC', '<i8', ()),
             ('MAXITERC', '<i8', ()), ('MINITERC', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHAR(object):
    name = 'CHAR'
    path = '/NASTRAN/INPUT/PARAMETER/PVT'
    dtype = [('NAME', 'S8', ()), ('VALUE', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CPLX(object):
    name = 'CPLX'
    path = '/NASTRAN/INPUT/PARAMETER/PVT'
    dtype = [('NAME', 'S8', ()), ('RE', '<f8', ()), ('IM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DOUBLE(object):
    name = 'DOUBLE'
    path = '/NASTRAN/INPUT/PARAMETER/PVT'
    dtype = [('NAME', 'S8', ()), ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INT(object):
    name = 'INT'
    path = '/NASTRAN/INPUT/PARAMETER/PVT'
    dtype = [('NAME', 'S8', ()), ('VALUE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPRELAX(object):
    name = 'SPRELAX'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID1', '<i8', ()), ('SID2', '<i8', ()), ('LIST2', '<i8', ()), ('DREF', '<f8', ()), ('LIST1', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SWLDPRM(object):
    name = 'SWLDPRM'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('TYPE', 'S8', ()), ('GSPROJ', '<f8', ()), ('PROJTOL', '<f8', ()), ('GSMOVE', '<i8', ()),
             ('CHKRUN', '<i8', ()), ('PRTSW', '<i8', ()), ('SAVSW', '<i8', ()), ('NREDIA', '<i8', ()),
             ('GSTOL', '<f8', ()), ('GMCHK', '<i8', ()), ('CNRAGLI', '<f8', ()), ('CNRAGLO', '<f8', ()),
             ('CSVOUT', '<i8', ()), ('MOVGAB', '<i8', ()), ('WMASS', '<i8', ()), ('SKIN', '<i8', ()),
             ('SCLSKIN', '<f8', ()), ('DRATIO', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRMCPL(object):
    name = 'TRMCPL'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('TID', '<i8', ()), ('CTYPE', 'S8', ()), ('PLTOL', '<f8', ()), ('GAPTOL1', '<f8', ()),
             ('GAPTOL2', '<f8', ()), ('GAPTOL3', '<f8', ()), ('GAPTOL4', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class T(object):
    name = 'T'
    path = '/NASTRAN/INPUT/PARAMETER/TSTEP'
    dtype = [('N', '<i8', ()), ('DT', '<f8', ()), ('NO', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARAMETER/TSTEP'
    dtype = [('SID', '<i8', ()), ('T_POS', '<i8', ()), ('T_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARAMETER/TSTEP/T']


@register_table
class TSTEPNL(object):
    name = 'TSTEPNL'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('SID', '<i8', ()), ('NDT', '<i8', ()), ('DT', '<f8', ()), ('NO', '<i8', ()), ('METHOD', '<i8', ()),
             ('KSTEP', '<i8', ()), ('MAXITER', '<i8', ()), ('CONV', '<i8', ()), ('EPSU', '<f8', ()),
             ('EPSP', '<f8', ()), ('EPSW', '<f8', ()), ('MAXDIV', '<i8', ()), ('MAXQN', '<i8', ()),
             ('MAXLS', '<i8', ()), ('FSTRESS', '<f8', ()), ('MAXBIS', '<i8', ()), ('ADJUST', '<i8', ()),
             ('MSTEP', '<i8', ()), ('RB', '<f8', ()), ('MAXR', '<f8', ()), ('UTOL', '<f8', ()), ('RTOLB', '<f8', ()),
             ('MINITER', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class UDNAME(object):
    name = 'UDNAME'
    path = '/NASTRAN/INPUT/PARAMETER'
    dtype = [('UID', '<i8', ()), ('NAME', 'S256', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACTIVAT(object):
    name = 'ACTIVAT'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('IPOST', '<i8', ()), ('ISET', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIS(object):
    name = 'EIS'
    path = '/NASTRAN/INPUT/PARTITION/AELIST'
    dtype = [('E', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/AELIST'
    dtype = [('SID', '<i8', ()), ('EIS_POS', '<i8', ()), ('EIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/AELIST/EIS']


@register_table
class CIS(object):
    name = 'CIS'
    path = '/NASTRAN/INPUT/PARTITION/AELISTC'
    dtype = [('CI', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/AELISTC'
    dtype = [('SID', '<i8', ()), ('CIS_POS', '<i8', ()), ('CIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/AELISTC/CIS']


@register_table
class AESURF(object):
    name = 'AESURF'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('CID1', '<i8', ()), ('ALID1', '<i8', ()), ('CID2', '<i8', ()),
             ('ALID2', '<i8', ()), ('EFF', '<f8', ()), ('LDW', '<i8', ()), ('CREFC', '<f8', ()), ('CREFS', '<f8', ()),
             ('S_PLLIM', '<i8', ()), ('PLLIM', '<f8', ()), ('S_PULIM', '<i8', ()), ('PULIM', '<f8', ()),
             ('S_HMLLIM', '<i8', ()), ('HMLLIM', '<f8', ()), ('S_HMULIM', '<i8', ()), ('HMULIM', '<f8', ()),
             ('S_TQLLIM', '<i8', ()), ('TQLLIM', '<i8', ()), ('S_TQULIM', '<i8', ()), ('TQULIM', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AESURFS(object):
    name = 'AESURFS'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('LABEL', 'S8', ()), ('LIST1', '<i8', ()), ('LIST2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ASET(object):
    name = 'ASET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/ASET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/ASET1'
    dtype = [('C', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/ASET1/IDS']


@register_table
class IDFS(object):
    name = 'IDFS'
    path = '/NASTRAN/INPUT/PARTITION/BDYLIST'
    dtype = [('IDF', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/BDYLIST'
    dtype = [('RHO', '<f8', ()), ('IDFS_POS', '<i8', ()), ('IDFS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/BDYLIST/IDFS']


@register_table
class GPS(object):
    name = 'GPS'
    path = '/NASTRAN/INPUT/PARTITION/BNDGRID'
    dtype = [('GP', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/BNDGRID'
    dtype = [('C', '<i8', ()), ('GPS_POS', '<i8', ()), ('GPS_LEN', '<i8', ()), ('GP1', '<i8', ()), ('GP2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/BNDGRID/GPS']


@register_table
class BSET(object):
    name = 'BSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/BSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/BSET1'
    dtype = [('C', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/BSET1/IDS']


@register_table
class GPS(object):
    name = 'GPS'
    path = '/NASTRAN/INPUT/PARTITION/CSUPER'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/CSUPER'
    dtype = [('SSID', '<i8', ()), ('PSID', '<i8', ()), ('GPS_POS', '<i8', ()), ('GPS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/CSUPER/GPS']


@register_table
class GPS(object):
    name = 'GPS'
    path = '/NASTRAN/INPUT/PARTITION/CSUPEXT'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/CSUPEXT'
    dtype = [('SEID', '<i8', ()), ('GPS_POS', '<i8', ()), ('GPS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/CSUPEXT/GPS']


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/PARTITION/CYJOIN'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/CYJOIN'
    dtype = [('SIDE', '<i8', ()), ('C', 'S8', ()), ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/CYJOIN/GRIDS']


@register_table
class CYSYM(object):
    name = 'CYSYM'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('NSEG', '<i8', ()), ('STYPE', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DEACTEL(object):
    name = 'DEACTEL'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('STRESS', '<i8', ()), ('STRAIN', '<i8', ()), ('IPOST', '<i8', ()), ('ISET', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIDS(object):
    name = 'EIDS'
    path = '/NASTRAN/INPUT/PARTITION/ELIST'
    dtype = [('EI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/ELIST'
    dtype = [('LID', '<i8', ()), ('EIDS_POS', '<i8', ()), ('EIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/ELIST/EIDS']


@register_table
class PANELS(object):
    name = 'PANELS'
    path = '/NASTRAN/INPUT/PARTITION/ERPPNL'
    dtype = [('NAME', 'S8', ()), ('SETID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/ERPPNL'
    dtype = [('PANELS_POS', '<i8', ()), ('PANELS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/ERPPNL/PANELS']


@register_table
class GIDS(object):
    name = 'GIDS'
    path = '/NASTRAN/INPUT/PARTITION/EXTRN'
    dtype = [('GID', '<i8', ()), ('C', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/EXTRN'
    dtype = [('GIDS_POS', '<i8', ()), ('GIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/EXTRN/GIDS']


@register_table
class FREEPT(object):
    name = 'FREEPT'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('IDF', '<i8', ()), ('IDP', '<i8', ()), ('PHI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDFS(object):
    name = 'IDFS'
    path = '/NASTRAN/INPUT/PARTITION/FSLIST'
    dtype = [('IDF', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/FSLIST'
    dtype = [('RHO', '<f8', ()), ('IDFS_POS', '<i8', ()), ('IDFS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/FSLIST/IDFS']


@register_table
class GRIDA(object):
    name = 'GRIDA'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('PGRID', '<i8', ()), ('PARTNAME', 'S64', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class OMIT(object):
    name = 'OMIT'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/OMIT1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/OMIT1'
    dtype = [('C', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/OMIT1/IDS']


@register_table
class OMITAX(object):
    name = 'OMITAX'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('RID', '<i8', ()), ('HID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PANELS(object):
    name = 'PANELS'
    path = '/NASTRAN/INPUT/PARTITION/PANEL'
    dtype = [('NAME', 'S8', ()), ('SETID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/PANEL'
    dtype = [('PANELS_POS', '<i8', ()), ('PANELS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/PANEL/PANELS']


@register_table
class QSET(object):
    name = 'QSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/QSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/QSET1'
    dtype = [('C', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/QSET1/IDS']


@register_table
class SETS(object):
    name = 'SETS'
    path = '/NASTRAN/INPUT/PARTITION/RADCAV'
    dtype = [('SET', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/RADCAV'
    dtype = [('ICAVITY', '<i8', ()), ('ELEAMB', '<i8', ()), ('SHADOW', '<i8', ()), ('SCALE', '<f8', ()),
             ('PRTPCH', '<i8', ()), ('NFECI', '<i8', ()), ('RMAX', '<f8', ()), ('NCOMP', '<i8', ()),
             ('SETS_POS', '<i8', ()), ('SETS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/RADCAV/SETS']


@register_table
class EIDS(object):
    name = 'EIDS'
    path = '/NASTRAN/INPUT/PARTITION/RADLST'
    dtype = [('EID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/RADLST'
    dtype = [('ICAVITY', '<i8', ()), ('MTXTYP', '<i8', ()), ('EIDS_POS', '<i8', ()), ('EIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/RADLST/EIDS']


@register_table
class MATRIX(object):
    name = 'MATRIX'
    path = '/NASTRAN/INPUT/PARTITION/RADMTX'
    dtype = [('F', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/RADMTX'
    dtype = [('ICAVITY', '<i8', ()), ('INDEX', '<i8', ()), ('MATRIX_POS', '<i8', ()), ('MATRIX_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/RADMTX/MATRIX']


@register_table
class RADSET(object):
    name = 'RADSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ICAVITY', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/RELEASE'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/RELEASE'
    dtype = [('SEID', '<i8', ()), ('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/RELEASE/IDS']


@register_table
class RVDOF(object):
    name = 'RVDOF'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/RVDOF1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/RVDOF1'
    dtype = [('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/RVDOF1/IDS']


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/PARTITION/SEBNDRY'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEBNDRY'
    dtype = [('SEIDA', '<i8', ()), ('SEIDB', '<i8', ()), ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEBNDRY/GRIDS']


@register_table
class SEBSET(object):
    name = 'SEBSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SEBSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEBSET1'
    dtype = [('SEID', '<i8', ()), ('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEBSET1/IDS']


@register_table
class SEBULK(object):
    name = 'SEBULK'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('TYPE', '<i8', ()), ('RSEID', '<i8', ()), ('METHOD', '<i8', ()), ('TOL', '<f8', ()),
             ('LOC', '<i8', ()), ('MEDIA', '<i8', ()), ('UNIT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/PARTITION/SECONCT'
    dtype = [('GA', '<i8', ()), ('GB', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SECONCT'
    dtype = [('SEIDA', '<i8', ()), ('SEIDB', '<i8', ()), ('TOL', '<f8', ()), ('LOC', '<i8', ()),
             ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SECONCT/GRIDS']


@register_table
class SECSET(object):
    name = 'SECSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SECSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SECSET1'
    dtype = [('SEID', '<i8', ()), ('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SECSET1/IDS']


@register_table
class EIDS(object):
    name = 'EIDS'
    path = '/NASTRAN/INPUT/PARTITION/SEELT'
    dtype = [('EID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEELT'
    dtype = [('SEID', '<i8', ()), ('EIDS_POS', '<i8', ()), ('EIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEELT/EIDS']


@register_table
class GRIDS(object):
    name = 'GRIDS'
    path = '/NASTRAN/INPUT/PARTITION/SEEXCLD'
    dtype = [('GA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEEXCLD'
    dtype = [('SEIDA', '<i8', ()), ('SEIDB', '<i8', ()), ('GRIDS_POS', '<i8', ()), ('GRIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEEXCLD/GRIDS']


@register_table
class SELABEL(object):
    name = 'SELABEL'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('LABEL', 'S56', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SELOC(object):
    name = 'SELOC'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('GA1', '<i8', ()), ('GA2', '<i8', ()), ('GA3', '<i8', ()), ('GB1', '<i8', ()),
             ('GB2', '<i8', ()), ('GB3', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SEMPLN(object):
    name = 'SEMPLN'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('FLAG', '<i8', ()), ('G1', '<i8', ()), ('G2', '<i8', ()), ('G3', '<i8', ()),
             ('CID', '<i8', ()), ('N1', '<f8', ()), ('N2', '<f8', ()), ('N3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SENQSET(object):
    name = 'SENQSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('NQSET', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GPS(object):
    name = 'GPS'
    path = '/NASTRAN/INPUT/PARTITION/SEQSEP'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEQSEP'
    dtype = [('SSID', '<i8', ()), ('PSID', '<i8', ()), ('GPS_POS', '<i8', ()), ('GPS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEQSEP/GPS']


@register_table
class SEQSET(object):
    name = 'SEQSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SEQSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEQSET1'
    dtype = [('SEID', '<i8', ()), ('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEQSET1/IDS']


@register_table
class GIDS(object):
    name = 'GIDS'
    path = '/NASTRAN/INPUT/PARTITION/SESET'
    dtype = [('G', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SESET'
    dtype = [('SEID', '<i8', ()), ('GIDS_POS', '<i8', ()), ('GIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SESET/GIDS']


@register_table
class GIDS(object):
    name = 'GIDS'
    path = '/NASTRAN/INPUT/PARTITION/SET1'
    dtype = [('G1', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SET1'
    dtype = [('SID', '<i8', ()), ('GIDS_POS', '<i8', ()), ('GIDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SET1/GIDS']


@register_table
class SET2(object):
    name = 'SET2'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SID', '<i8', ()), ('MACRO', '<i8', ()), ('SP1', '<f8', ()), ('SP2', '<f8', ()), ('CH1', '<f8', ()),
             ('CH2', '<f8', ()), ('ZMAX', '<f8', ()), ('ZMIN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SET3'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SET3'
    dtype = [('SID', '<i8', ()), ('DES', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SET3/IDS']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SET4'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SET4'
    dtype = [('SID', '<i8', ()), ('CLASS', 'S8', ()), ('TYPE', 'S8', ()), ('IDS_POS', '<i8', ()),
             ('IDS_LEN', '<i8', ()), ('FLAG_ALL', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SET4/IDS']


@register_table
class SEUPIS(object):
    name = 'SEUPIS'
    path = '/NASTRAN/INPUT/PARTITION/SETREE'
    dtype = [('SEUPI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SETREE'
    dtype = [('SEID', '<i8', ()), ('SEUPIS_POS', '<i8', ()), ('SEUPIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SETREE/SEUPIS']


@register_table
class SEUSET(object):
    name = 'SEUSET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SEID', '<i8', ()), ('SNAME', '<i8', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SEUSET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SEUSET1'
    dtype = [('SEID', '<i8', ()), ('SNAME', '<i8', ()), ('C', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ()),
             ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SEUSET1/IDS']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/SLBDY'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/SLBDY'
    dtype = [('RHO', '<f8', ()), ('M', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/SLBDY/IDS']


@register_table
class USET(object):
    name = 'USET'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('SNAME', 'S4', ()), ('ID', '<i8', ()), ('C', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PARTITION/USET1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/USET1'
    dtype = [('SNAME', 'S4', ()), ('C', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('ID1', '<i8', ()),
             ('ID2', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/USET1/IDS']


@register_table
class VIEW(object):
    name = 'VIEW'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('IVIEW', '<i8', ()), ('ICAVITY', '<i8', ()), ('SHADE', '<i8', ()), ('NB', '<i8', ()), ('NG', '<i8', ()),
             ('DISLIN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VIEW3D(object):
    name = 'VIEW3D'
    path = '/NASTRAN/INPUT/PARTITION'
    dtype = [('ICAVITY', '<i8', ()), ('GITB', '<i8', ()), ('GIPS', '<i8', ()), ('CIER', '<i8', ()), ('ETOL', '<f8', ()),
             ('ZTOL', '<f8', ()), ('WTOL', '<f8', ()), ('RADCHK', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WEIDS(object):
    name = 'WEIDS'
    path = '/NASTRAN/INPUT/PARTITION/WETSURF'
    dtype = [('WEID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PARTITION/WETSURF'
    dtype = [('WSID', '<i8', ()), ('WTAG', 'S8', ()), ('WEIDS_POS', '<i8', ()), ('WEIDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PARTITION/WETSURF/WEIDS']


@register_table
class MFLUID(object):
    name = 'MFLUID'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('SID', '<i8', ()), ('CID', '<i8', ()), ('ZFR', '<f8', ()), ('RHO', '<f8', ()), ('ELIST1', '<i8', ()),
             ('ELIST2', '<i8', ()), ('PLANE1', '<i8', ()), ('PLANE2', '<i8', ()), ('RMAX', '<f8', ()),
             ('FMEXACT', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PROPERTY/NSM'
    dtype = [('ID', '<i8', ()), ('VALUE', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/NSM'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S8', ()), ('PROP', 'S4', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/NSM/IDS']


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/PROPERTY/NSM1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/PROPERTY/NSM1'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU_BY(object):
    name = 'THRU_BY'
    path = '/NASTRAN/INPUT/PROPERTY/NSM1'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('N', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/NSM1'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S8', ()), ('PROP', 'S8', ()), ('VALUE', '<f8', ()), ('ALL', '<i8', ()),
             ('LIST_POS', '<i8', ()), ('LIST_LEN', '<i8', ()), ('THRU_POS', '<i8', ()), ('THRU_LEN', '<i8', ()),
             ('THRUBY_POS', '<i8', ()), ('THRUBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/NSM1/IDLIST', '/NASTRAN/INPUT/PROPERTY/NSM1/THRU',
                 '/NASTRAN/INPUT/PROPERTY/NSM1/THRU_BY']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PROPERTY/NSMADD'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/NSMADD'
    dtype = [('SID', '<i8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/NSMADD/IDS']


@register_table
class IDS(object):
    name = 'IDS'
    path = '/NASTRAN/INPUT/PROPERTY/NSML'
    dtype = [('ID', '<i8', ()), ('VALUE', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/NSML'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S8', ()), ('IDS_POS', '<i8', ()), ('IDS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/NSML/IDS']


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/PROPERTY/NSML1'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/PROPERTY/NSML1'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU_BY(object):
    name = 'THRU_BY'
    path = '/NASTRAN/INPUT/PROPERTY/NSML1'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('N', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/NSML1'
    dtype = [('SID', '<i8', ()), ('TYPE', 'S8', ()), ('VALUE', '<f8', ()), ('ALL', '<i8', ()), ('LIST_POS', '<i8', ()),
             ('LIST_LEN', '<i8', ()), ('THRU_POS', '<i8', ()), ('THRU_LEN', '<i8', ()), ('THRUBY_POS', '<i8', ()),
             ('THRUBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/NSML1/IDLIST', '/NASTRAN/INPUT/PROPERTY/NSML1/THRU',
                 '/NASTRAN/INPUT/PROPERTY/NSML1/THRU_BY']


@register_table
class PAABSF(object):
    name = 'PAABSF'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('TZREID', '<i8', ()), ('TZMID', '<i8', ()), ('S', '<f8', ()), ('A', '<f8', ()),
             ('B', '<f8', ()), ('K', '<f8', ()), ('RHOC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PACABS(object):
    name = 'PACABS'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('SYNTH', '<i8', ()), ('TID1', '<i8', ()), ('TID2', '<i8', ()), ('TID3', '<i8', ()),
             ('TESTAR', '<f8', ()), ('CUTFR', '<f8', ()), ('B', '<f8', ()), ('K', '<f8', ()), ('M', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PACBAR(object):
    name = 'PACBAR'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MBACK', '<f8', ()), ('MSEPTM', '<f8', ()), ('FRESON', '<f8', ()),
             ('KRESON', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PACINF(object):
    name = 'PACINF'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('RIO', '<i8', ()), ('XP', '<f8', ()), ('YP', '<f8', ()),
             ('ZP', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PAERO1(object):
    name = 'PAERO1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('B1', '<i8', ()), ('B2', '<i8', ()), ('B3', '<i8', ()), ('B4', '<i8', ()),
             ('B5', '<i8', ()), ('B6', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PAERO2(object):
    name = 'PAERO2'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('ORIENT', 'S4', ()), ('WIDTH', '<f8', ()), ('AR', '<f8', ()), ('LRSB', '<i8', ()),
             ('LRIB', '<i8', ()), ('LTH1', '<i8', ()), ('LTH2', '<i8', ()), ('THI1', '<i8', ()), ('THN1', '<i8', ()),
             ('THI2', '<i8', ()), ('THN2', '<i8', ()), ('THI3', '<i8', ()), ('THN3', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PAERO3(object):
    name = 'PAERO3'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('NBOX', '<i8', ()), ('NCTRL', '<i8', ()), ('X5', '<f8', ()), ('Y5', '<f8', ()),
             ('X6', '<f8', ()), ('Y6', '<f8', ()), ('X7', '<f8', ()), ('Y7', '<f8', ()), ('X8', '<f8', ()),
             ('Y8', '<f8', ()), ('X9', '<f8', ()), ('Y9', '<f8', ()), ('X10', '<f8', ()), ('Y10', '<f8', ()),
             ('X11', '<f8', ()), ('Y11', '<f8', ()), ('X12', '<f8', ()), ('Y12', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FIELDS(object):
    name = 'FIELDS'
    path = '/NASTRAN/INPUT/PROPERTY/PAERO4'
    dtype = [('DOCI', '<f8', ()), ('CAOCI', '<f8', ()), ('GAPOCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PAERO4'
    dtype = [('PID', '<i8', ()), ('CLA', '<i8', ()), ('LCLA', '<i8', ()), ('CIRC', '<i8', ()), ('LCIRC', '<i8', ()),
             ('FIELDS_POS', '<i8', ()), ('FIELDS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PAERO4/FIELDS']


@register_table
class CAOCIS(object):
    name = 'CAOCIS'
    path = '/NASTRAN/INPUT/PROPERTY/PAERO5'
    dtype = [('CAOCI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PAERO5'
    dtype = [('PID', '<i8', ()), ('NALPHA', '<i8', ()), ('LALPHA', '<i8', ()), ('NXIS', '<i8', ()), ('LXIS', '<i8', ()),
             ('NTAUS', '<i8', ()), ('LTAUS', '<i8', ()), ('CAOCIS_POS', '<i8', ()), ('CAOCIS_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PAERO5/CAOCIS']


@register_table
class PAXISYM(object):
    name = 'PAXISYM'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('INT', 'S4', ()), ('T1', '<f8', ()), ('T2', '<f8', ()),
             ('ANAL', 'S8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PAXSYMH(object):
    name = 'PAXSYMH'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('CID', '<i8', ()), ('NHARM', '<i8', ()), ('INT', '<i8', ()),
             ('ROTORID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PBAR(object):
    name = 'PBAR'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('A', '<f8', ()), ('I1', '<f8', ()), ('I2', '<f8', ()),
             ('J', '<f8', ()), ('NSM', '<f8', ()), ('FE', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()),
             ('D1', '<f8', ()), ('D2', '<f8', ()), ('E1', '<f8', ()), ('E2', '<f8', ()), ('F1', '<f8', ()),
             ('F2', '<f8', ()), ('K1', '<f8', ()), ('K2', '<f8', ()), ('I12', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INFO(object):
    name = 'INFO'
    path = '/NASTRAN/INPUT/PROPERTY/PBARL'
    dtype = [('VALUE', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBARL'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('GROUP', 'S8', ()), ('TYPE', 'S8', ()), ('INFO_POS', '<i8', ()),
             ('INFO_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBARL/INFO']


@register_table
class PBARN1(object):
    name = 'PBARN1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('RADIUS', '<f8', ()), ('SECT', 'S4', ()), ('IANAL', '<i8', ()),
             ('BEH2', 'S4', ()), ('INT2', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBCOMP'
    dtype = [('Y', '<f8', ()), ('Z', '<f8', ()), ('C', '<f8', ()), ('MID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBCOMP'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('A', '<f8', ()), ('I1', '<f8', ()), ('I2', '<f8', ()),
             ('I12', '<f8', ()), ('J', '<f8', ()), ('NSM', '<f8', ()), ('K1', '<f8', ()), ('K2', '<f8', ()),
             ('M1', '<f8', ()), ('M2', '<f8', ()), ('N1', '<f8', ()), ('N2', '<f8', ()), ('NSECT', '<i8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBCOMP/SECTION']


@register_table
class PBEAM(object):
    name = 'PBEAM'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('NSEGS', '<i8', ()), ('CCF', '<i8', ()), ('CWELD', '<i8', ()), (
    'SECTION',
    [('SO', '<f8', ()), ('XXB', '<f8', ()), ('A', '<f8', ()), ('I1', '<f8', ()), ('I2', '<f8', ()), ('I12', '<f8', ()),
     ('J', '<f8', ()), ('NSM', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()), ('D1', '<f8', ()), ('D2', '<f8', ()),
     ('E1', '<f8', ()), ('E2', '<f8', ()), ('F1', '<f8', ()), ('F2', '<f8', ())], (11,)), ('K1', '<f8', ()),
             ('K2', '<f8', ()), ('S1', '<f8', ()), ('S2', '<f8', ()), ('NSIA', '<f8', ()), ('NSIB', '<f8', ()),
             ('CWA', '<f8', ()), ('CWB', '<f8', ()), ('M1A', '<f8', ()), ('M2A', '<f8', ()), ('M1B', '<f8', ()),
             ('M2B', '<f8', ()), ('N1A', '<f8', ()), ('N2A', '<f8', ()), ('N1B', '<f8', ()), ('N2B', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PBEAM3(object):
    name = 'PBEAM3'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('AA', '<f8', ()), ('IZA', '<f8', ()), ('IYA', '<f8', ()),
             ('IYZA', '<f8', ()), ('JA', '<f8', ()), ('NSMA', '<f8', ()), ('CYA', '<f8', ()), ('CZA', '<f8', ()),
             ('DYA', '<f8', ()), ('DZA', '<f8', ()), ('EYA', '<f8', ()), ('EZA', '<f8', ()), ('FYA', '<f8', ()),
             ('FZA', '<f8', ()), ('SECTION', [('SO', '<i8', ()), ('A', '<f8', ()), ('IZ', '<f8', ()), ('IY', '<f8', ()),
                                              ('IYZ', '<f8', ()), ('J', '<f8', ()), ('NSM', '<f8', ()),
                                              ('CY', '<f8', ()), ('CZ', '<f8', ()), ('DY', '<f8', ()),
                                              ('DZ', '<f8', ()), ('EY', '<f8', ()), ('EZ', '<f8', ()),
                                              ('FY', '<f8', ()), ('FZ', '<f8', ())], (2,)), ('KY', '<f8', ()),
             ('KZ', '<f8', ()), ('NYA', '<f8', ()), ('NZA', '<f8', ()), ('NYB', '<f8', ()), ('NZB', '<f8', ()),
             ('NYC', '<f8', ()), ('NZC', '<f8', ()), ('MYA', '<f8', ()), ('MZA', '<f8', ()), ('MYB', '<f8', ()),
             ('MZB', '<f8', ()), ('MYC', '<f8', ()), ('MZC', '<f8', ()), ('NSIA', '<f8', (3,)), ('NSIB', '<f8', (3,)),
             ('NSIC', '<f8', (3,)), ('CW', '<f8', (3,)), ('STRESS', '<i8', ()), ('WIJ', '<f8', (36,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DIMS(object):
    name = 'DIMS'
    path = '/NASTRAN/INPUT/PROPERTY/PBEAML'
    dtype = [('DIM', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBEAML'
    dtype = [('SO', '<f8', ()), ('RDIST', '<f8', ()), ('DIMS_POS', '<i8', ()), ('DIMS_LEN', '<i8', ()),
             ('NSM', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBEAML'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('GROUP', 'S8', ()), ('TYPE', 'S8', ()),
             ('SECTION_POS', '<i8', ()), ('SECTION_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBEAML/DIMS', '/NASTRAN/INPUT/PROPERTY/PBEAML/SECTION']


@register_table
class PBEMN1(object):
    name = 'PBEMN1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('RADIUS', '<f8', ()), ('SECT', 'S4', ()), ('IANAL', '<i8', ()),
             ('BEH2', 'S4', ()), ('INT2', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PBEND(object):
    name = 'PBEND'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('A', '<f8', ()), ('I1', '<f8', ()), ('I2', '<f8', ()),
             ('J', '<f8', ()), ('FSI', '<i8', ()), ('RM', '<f8', ()), ('T', '<f8', ()), ('P', '<f8', ()),
             ('RB', '<f8', ()), ('THETAB', '<f8', ()), ('C1', '<f8', ()), ('C2', '<f8', ()), ('D1', '<f8', ()),
             ('D2', '<f8', ()), ('E1', '<f8', ()), ('E2', '<f8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()),
             ('K1', '<f8', ()), ('K2', '<f8', ()), ('NSM', '<f8', ()), ('RC', '<f8', ()), ('ZC', '<f8', ()),
             ('DELTAN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ING(object):
    name = 'ING'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/GS'
    dtype = [('EXTID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INP(object):
    name = 'INP'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/GS'
    dtype = [('EXTID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/GS'
    dtype = [('STATID', '<i8', ()), ('SRECL', '<i8', ()), ('SO', '<f8', ()), ('XXB', '<f8', ()), ('NSM', '<f8', ()),
             ('OUTP', '<i8', ()), ('OUTG', '<i8', ()), ('OUTM', '<i8', ()), ('OHEIGHT', '<f8', ()),
             ('OWIDTH', '<f8', ()), ('INP_POS', '<i8', ()), ('INP_LEN', '<i8', ()), ('ING_POS', '<i8', ()),
             ('ING_LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BRG(object):
    name = 'BRG'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTBRGID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BRP(object):
    name = 'BRP'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTBRPID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORE(object):
    name = 'CORE'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTCORID', '<i8', ()), ('EXTCPCID', '<i8', ()), ('TYPE', '<i8', ()), ('ID1', '<i8', ()),
             ('ID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEIGHT(object):
    name = 'HEIGHT'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTHGTID', '<i8', ()), ('HGPTID1', '<i8', ()), ('HGPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LAYER(object):
    name = 'LAYER'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTLAYID', '<i8', ()), ('EXTLPCID', '<i8', ()), ('EXTLSTID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LENGTH(object):
    name = 'LENGTH'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTLENID', '<i8', ()), ('LNPTID1', '<i8', ()), ('LNPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('STATID', '<i8', ()), ('SRECL', '<i8', ()), ('SO', '<f8', ()), ('XXB', '<f8', ()), ('NSM', '<f8', ()),
             ('OUTP', '<i8', ()), ('OUTG', '<i8', ()), ('OUTM', '<i8', ()), ('OHEIGHT', '<f8', ()),
             ('OWIDTH', '<f8', ()), ('BRP_POS', '<i8', ()), ('BRP_LEN', '<i8', ()), ('BRG_POS', '<i8', ()),
             ('BRG_LEN', '<i8', ()), ('HEIGHT_POS', '<i8', ()), ('HEIGHT_LEN', '<i8', ()), ('WIDTH_POS', '<i8', ()),
             ('WIDTH_LEN', '<i8', ()), ('LENGTH_POS', '<i8', ()), ('LENGTH_LEN', '<i8', ()),
             ('THICKNESS_POS', '<i8', ()), ('THICKNESS_LEN', '<i8', ()), ('CORE_POS', '<i8', ()),
             ('CORE_LEN', '<i8', ()), ('LAYER_POS', '<i8', ()), ('LAYER_LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class THICKNESS(object):
    name = 'THICKNESS'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTTHKID', '<i8', ()), ('THICK', '<f8', ()), ('TYPE', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WIDTH(object):
    name = 'WIDTH'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP'
    dtype = [('EXTWIDID', '<i8', ()), ('WDPTID1', '<i8', ()), ('WDPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBMSECT'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('FORM', '<i8', ()), ('SECTION_POS', '<i8', ()),
             ('SECTION_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBMSECT/GS', '/NASTRAN/INPUT/PROPERTY/PBMSECT/OP_CP']


@register_table
class ING(object):
    name = 'ING'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/GS'
    dtype = [('EXTID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INP(object):
    name = 'INP'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/GS'
    dtype = [('EXTID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/GS'
    dtype = [('OHEIGHT', '<f8', ()), ('OWIDTH', '<f8', ()), ('INP_POS', '<i8', ()), ('INP_LEN', '<i8', ()),
             ('ING_POS', '<i8', ()), ('ING_LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BRG(object):
    name = 'BRG'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTBRGID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BRP(object):
    name = 'BRP'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTBRPID', '<i8', ()), ('ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CORE(object):
    name = 'CORE'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTCORID', '<i8', ()), ('EXTCPCID', '<i8', ()), ('TYPE', '<i8', ()), ('ID1', '<i8', ()),
             ('ID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEIGHT(object):
    name = 'HEIGHT'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTHGTID', '<i8', ()), ('HGPTID1', '<i8', ()), ('HGPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LAYER(object):
    name = 'LAYER'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTLAYID', '<i8', ()), ('EXTLPCID', '<i8', ()), ('EXTLSTID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LENGTH(object):
    name = 'LENGTH'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTLENID', '<i8', ()), ('LNPTID1', '<i8', ()), ('LNPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SECTION(object):
    name = 'SECTION'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('OHEIGHT', '<f8', ()), ('OWIDTH', '<f8', ()), ('BRP_POS', '<i8', ()), ('BRP_LEN', '<i8', ()),
             ('BRG_POS', '<i8', ()), ('BRG_LEN', '<i8', ()), ('HEIGHT_POS', '<i8', ()), ('HEIGHT_LEN', '<i8', ()),
             ('WIDTH_POS', '<i8', ()), ('WIDTH_LEN', '<i8', ()), ('LENGTH_POS', '<i8', ()), ('LENGTH_LEN', '<i8', ()),
             ('THICKNESS_POS', '<i8', ()), ('THICKNESS_LEN', '<i8', ()), ('CORE_POS', '<i8', ()),
             ('CORE_LEN', '<i8', ()), ('LAYER_POS', '<i8', ()), ('LAYER_LEN', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class THICKNESS(object):
    name = 'THICKNESS'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTTHKID', '<i8', ()), ('THICK', '<f8', ()), ('TYPE', '<i8', ()), ('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WIDTH(object):
    name = 'WIDTH'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP'
    dtype = [('EXTWIDID', '<i8', ()), ('WDPTID1', '<i8', ()), ('WDPTID2', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBRSECT'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('FORM', '<i8', ()), ('RECLEN', '<i8', ()), ('NSM', '<f8', ()),
             ('OUTP', '<i8', ()), ('OUTG', '<i8', ()), ('OUTM', '<i8', ()), ('SECTION_POS', '<i8', ()),
             ('SECTION_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBRSECT/GS', '/NASTRAN/INPUT/PROPERTY/PBRSECT/OP_CP']


@register_table
class PBUSH(object):
    name = 'PBUSH'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('K', '<f8', (6,)), ('B', '<f8', (6,)), ('GE', '<f8', (6,)), ('SA', '<f8', ()),
             ('ST', '<f8', ()), ('EA', '<f8', ()), ('ET', '<f8', ()), ('M', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PBUSH1D(object):
    name = 'PBUSH1D'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('K', '<f8', ()), ('C', '<f8', ()), ('M', '<f8', ()), ('ALPHA', '<f8', ()),
             ('SA', '<f8', ()), ('EA', '<f8', ()), ('TYPEA', '<i8', ()), ('CVT', '<f8', ()), ('CVC', '<f8', ()),
             ('EXPVT', '<f8', ()), ('EXPVC', '<f8', ()), ('IDTSU', '<i8', ()), ('IDTCU', '<i8', ()),
             ('IDTSUD', '<i8', ()), ('IDCSUD', '<i8', ()), ('TYPES', '<i8', ()), ('IDTS', '<i8', ()),
             ('IDCS', '<i8', ()), ('IDTDU1', '<i8', ()), ('IDCDU1', '<i8', ()), ('TYPED', '<i8', ()),
             ('IDTD1', '<i8', ()), ('IDTD2', '<i8', ()), ('IDTDV1', '<i8', ()), ('IDCDV1', '<i8', ()),
             ('TYPEG', '<i8', ()), ('IDTG', '<i8', ()), ('IDCG', '<i8', ()), ('IDTDU2', '<i8', ()),
             ('IDCDU2', '<i8', ()), ('IDTDV2', '<i8', ()), ('IDCDV2', '<i8', ()), ('TYPEF', '<i8', ()),
             ('IDTF', '<i8', ()), ('IDCF', '<i8', ()), ('UT', '<f8', ()), ('UC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CROSS(object):
    name = 'CROSS'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('K12', '<f8', ()), ('K21', '<f8', ()), ('B12', '<f8', ()), ('B21', '<f8', ()), ('M12', '<f8', ()),
             ('M21', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class FORCE(object):
    name = 'FORCE'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('FTBEQ', '<i8', ()), ('TIDF1', '<i8', ()), ('TIDF2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class RGAP(object):
    name = 'RGAP'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('TABK', '<i8', ()), ('TABB', '<i8', ()), ('TABG', '<i8', ()), ('TABU', '<i8', ()), ('RADIUS', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SDM(object):
    name = 'SDM'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('OPTTYP', '<i8', ()), ('TBEQ', '<i8', ()), ('TEID11', '<i8', ()), ('TEID22', '<i8', ()),
             ('TEID12', '<i8', ()), ('TEID21', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SQUEEZE(object):
    name = 'SQUEEZE'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('BDIA', '<f8', ()), ('BLEN', '<f8', ()), ('BCLR', '<f8', ()), ('SOLN', '<i8', ()), ('VISCO', '<f8', ()),
             ('PVAPCO', '<f8', ()), ('NPORT', '<i8', ()), ('PRES1', '<f8', ()), ('THETA1', '<f8', ()),
             ('PRES2', '<f8', ()), ('THETA2', '<f8', ()), ('OFFSET1', '<f8', ()), ('OFFSET2', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PBUSH2D'
    dtype = [('PID', '<i8', ()), ('K1', '<f8', ()), ('K2', '<f8', ()), ('B1', '<f8', ()), ('B2', '<f8', ()),
             ('M1', '<f8', ()), ('M2', '<f8', ()), ('FLAG', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PBUSH2D/CROSS', '/NASTRAN/INPUT/PROPERTY/PBUSH2D/FORCE',
                 '/NASTRAN/INPUT/PROPERTY/PBUSH2D/RGAP', '/NASTRAN/INPUT/PROPERTY/PBUSH2D/SDM',
                 '/NASTRAN/INPUT/PROPERTY/PBUSH2D/SQUEEZE']


@register_table
class PBUSHT(object):
    name = 'PBUSHT'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('TKID', '<i8', (6,)), ('TBID', '<i8', (6,)), ('TGEID', '<i8', (6,)),
             ('TKNID', '<i8', (6,)), ('FDC', 'S8', ()), ('FUSE', '<i8', ()), ('DIR', '<i8', ()), ('OPTION', 'S8', ()),
             ('LOWER', '<f8', ()), ('UPPER', '<f8', ()), ('FRATE', '<f8', ()), ('LRGR', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PCOHE(object):
    name = 'PCOHE'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('INT', '<i8', ()), ('T', '<f8', ()), ('COHEOUT', 'S4', ()),
             ('SECANT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLY(object):
    name = 'PLY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMP'
    dtype = [('MID', '<i8', ()), ('T', '<f8', ()), ('THETA', '<f8', ()), ('SOUT', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMP'
    dtype = [('PID', '<i8', ()), ('NPLIES', '<i8', ()), ('Z0', '<f8', ()), ('NSM', '<f8', ()), ('SB', '<f8', ()),
             ('FT', '<i8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('PLY_POS', '<i8', ()), ('PLY_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PCOMP/PLY']


@register_table
class IDLIST(object):
    name = 'IDLIST'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPF'
    dtype = [('ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRU(object):
    name = 'THRU'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPF'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THRUBY(object):
    name = 'THRUBY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPF'
    dtype = [('ID1', '<i8', ()), ('ID2', '<i8', ()), ('N', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPF'
    dtype = [('INT', '<i8', ()), ('ALL', '<i8', ()), ('LIST_POS', '<i8', ()), ('LIST_LEN', '<i8', ()),
             ('THRU_POS', '<i8', ()), ('THRU_LEN', '<i8', ()), ('THRUBY_POS', '<i8', ()), ('THRUBY_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PCOMPF/IDLIST', '/NASTRAN/INPUT/PROPERTY/PCOMPF/THRU',
                 '/NASTRAN/INPUT/PROPERTY/PCOMPF/THRUBY']


@register_table
class PLY(object):
    name = 'PLY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPG'
    dtype = [('GPLYID', '<i8', ()), ('MID', '<i8', ()), ('THICK', '<f8', ()), ('THETA', '<f8', ()), ('SOUT', '<i8', ()),
             ('MIDMTX', '<i8', ()), ('VF', '<f8', ()), ('VV', '<f8', ()), ('CTEMP', '<f8', ()), ('MOIST', '<f8', ()),
             ('CRIT', 'S4', ()), ('NFTI', '<i8', ()), ('FTI', 'S96', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPG'
    dtype = [('PID', '<i8', ()), ('NPLIES', '<i8', ()), ('Z0', '<f8', ()), ('NSM', '<f8', ()), ('SB', '<f8', ()),
             ('FT', '<i8', ()), ('TREF', '<f8', ()), ('GE', '<f8', ()), ('MICRO', 'S4', ()), ('PLY_POS', '<i8', ()),
             ('PLY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PCOMPG/PLY']


@register_table
class PLY(object):
    name = 'PLY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPLS'
    dtype = [('GPLYID', '<i8', ()), ('MID', '<i8', ()), ('THICK', '<f8', ()), ('THETA', '<f8', ()), ('SOUT', 'S4', ()),
             ('MIDMTX', '<i8', ()), ('VF', '<f8', ()), ('VV', '<f8', ()), ('CTEMP', '<f8', ()), ('MOIST', '<f8', ()),
             ('CRIT', 'S4', ()), ('NFTI', '<i8', ()), ('FTI', 'S96', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PCOMPLS'
    dtype = [('PID', '<i8', ()), ('DIRECT', '<i8', ()), ('CORDM', '<i8', ()), ('SB', '<f8', ()), ('IANAL', '<i8', ()),
             ('MICRO', 'S4', ()), ('BEH8', 'S4', ()), ('INT8', 'S4', ()), ('BEH8H', 'S4', ()), ('INT8H', 'S4', ()),
             ('BEH20', 'S4', ()), ('INT20', 'S4', ()), ('BEH20H', 'S4', ()), ('INT20H', 'S4', ()), ('NPLY', '<i8', ()),
             ('PLY_POS', '<i8', ()), ('PLY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PCOMPLS/PLY']


@register_table
class PCONEAX(object):
    name = 'PCONEAX'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID1', '<i8', ()), ('T1', '<f8', ()), ('MID2', '<i8', ()), ('I', '<f8', ()),
             ('MID3', '<i8', ()), ('T2', '<f8', ()), ('NSM', '<f8', ()), ('Z1', '<f8', ()), ('Z2', '<f8', ()),
             ('PHI', '<f8', (14,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PCONV(object):
    name = 'PCONV'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('FORM', '<i8', ()), ('EXPF', '<f8', ()), ('FTYPE', '<i8', ()),
             ('HCF1', '<f8', ()), ('HCF2', '<f8', ()), ('HCF3', '<f8', ()), ('HCF4', '<f8', ()), ('HCF5', '<f8', ()),
             ('HCF6', '<f8', ()), ('HCF7', '<f8', ()), ('HCF8', '<f8', ()), ('TID', '<i8', ()), ('CHLEN', '<f8', ()),
             ('GIDIN', '<i8', ()), ('CE', '<i8', ()), ('E', '<f8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PCONV1(object):
    name = 'PCONV1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PCONID', '<i8', ()), ('CORRID', '<i8', ()), ('MID', '<i8', ()), ('MDT', '<f8', ()), ('VEL', '<f8', ()),
             ('LEN', '<f8', ()), ('FLW', '<f8', ()), ('MDTTYPE', '<i8', ()), ('VELTYPE', '<i8', ()),
             ('LENTYPE', '<i8', ()), ('FLWTYPE', '<i8', ()), ('MDTT', '<i8', ()), ('VELT', '<i8', ()),
             ('LENT', '<i8', ()), ('FLWT', '<i8', ()), ('CNT', '<i8', ()), ('CI', '<f8', (24,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PCONVM(object):
    name = 'PCONVM'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('FORM', '<i8', ()), ('FLAG', '<i8', ()), ('COEF', '<f8', ()),
             ('EXPR', '<f8', ()), ('EXPPI', '<f8', ()), ('EXPPO', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PDAMP(object):
    name = 'PDAMP'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('B', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PDAMP5(object):
    name = 'PDAMP5'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('B', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PDAMPT(object):
    name = 'PDAMPT'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('TBID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PELAS(object):
    name = 'PELAS'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('K', '<f8', ()), ('GE', '<f8', ()), ('S', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PELAST(object):
    name = 'PELAST'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('TKID', '<i8', ()), ('TGEID', '<i8', ()), ('TKNID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PFAST(object):
    name = 'PFAST'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('D', '<f8', ()), ('CONNBEH', '<i8', ()), ('CONNTYPE', '<i8', ()),
             ('EXTCON', '<i8', ()), ('CONDTYPE', '<i8', ()), ('WELDTYPE', '<i8', ()), ('MINLEN', '<f8', ()),
             ('MAXLEN', '<f8', ()), ('GMCHK', '<i8', ()), ('SPCGS', '<i8', ()), ('CMASS', '<f8', ()), ('GE', '<f8', ()),
             ('MCID', '<i8', ()), ('MFLAG', '<i8', ()), ('KT', '<f8', (3,)), ('KR', '<f8', (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PGAP(object):
    name = 'PGAP'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('UO', '<f8', ()), ('FO', '<f8', ()), ('KA', '<f8', ()), ('KB', '<f8', ()),
             ('KT', '<f8', ()), ('MU1', '<f8', ()), ('MU2', '<f8', ()), ('TMAX', '<f8', ()), ('MAR', '<f8', ()),
             ('TRMIN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PHBDY(object):
    name = 'PHBDY'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('AF', '<f8', ()), ('D1', '<f8', ()), ('D2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLY(object):
    name = 'PLY'
    path = '/NASTRAN/INPUT/PROPERTY/PLCOMP'
    dtype = [('GPLYID', '<i8', ()), ('MID', '<i8', ()), ('THICK', '<f8', ()), ('THETA', '<f8', ()), ('SOUT', 'S4', ()),
             ('MIDMTX', '<i8', ()), ('VF', '<f8', ()), ('VV', '<f8', ()), ('CTEMP', '<f8', ()), ('MOIST', '<f8', ()),
             ('CRIT', 'S4', ()), ('NFTI', '<i8', ()), ('FTI', 'S96', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/PLCOMP'
    dtype = [('PID', '<i8', ()), ('DIRECT', '<i8', ()), ('THICKOP', '<f8', ()), ('SB', '<f8', ()), ('IANAL', '<i8', ()),
             ('MICRO', 'S4', ()), ('BEH4', 'S4', ()), ('INT4', 'S4', ()), ('BEH4H', 'S4', ()), ('INT4H', 'S4', ()),
             ('BEH8', 'S4', ()), ('INT8', 'S4', ()), ('BEH8H', 'S4', ()), ('INT8H', 'S4', ()), ('NPLY', '<i8', ()),
             ('PLY_POS', '<i8', ()), ('PLY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/PLCOMP/PLY']


@register_table
class PLPLANE(object):
    name = 'PLPLANE'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('CID', '<i8', ()), ('STR', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PLSOLID(object):
    name = 'PLSOLID'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('STR', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PMASS(object):
    name = 'PMASS'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('M', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PROD(object):
    name = 'PROD'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('A', '<f8', ()), ('J', '<f8', ()), ('C', '<f8', ()),
             ('NSM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRODN1(object):
    name = 'PRODN1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('IANAL', '<i8', ()), ('BEH2', 'S4', ()), ('INT2', 'S4', ()),
             ('BEH2H', 'S4', ()), ('INT2H', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSEAM(object):
    name = 'PSEAM'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('TYPE', '<i8', ()), ('WIDTH', '<f8', ()), ('THK', '<f8', ()),
             ('IN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSHEAR(object):
    name = 'PSHEAR'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('T', '<f8', ()), ('NSM', '<f8', ()), ('F1', '<f8', ()),
             ('F2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSHEARN(object):
    name = 'PSHEARN'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('IANAL', '<i8', ()), ('BEH4', 'S4', ()), ('INT4', 'S4', ()),
             ('BEH4H', 'S4', ()), ('INT4H', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSHELL(object):
    name = 'PSHELL'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID1', '<i8', ()), ('T', '<f8', ()), ('MID2', '<i8', ()), ('BK', '<f8', ()),
             ('MID3', '<i8', ()), ('TS', '<f8', ()), ('NSM', '<f8', ()), ('Z1', '<f8', ()), ('Z2', '<f8', ()),
             ('MID4', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSHLN1(object):
    name = 'PSHLN1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID1', '<i8', ()), ('MID2', '<i8', ()), ('SMEAR', 'S4', ()), ('IANAL', '<i8', ()),
             ('BEH3', 'S4', ()), ('INT3', 'S4', ()), ('BEH3H', 'S4', ()), ('INT3H', 'S4', ()), ('TMPDIST3', '<i8', ()),
             ('BEH4', 'S4', ()), ('INT4', 'S4', ()), ('BEH4H', 'S4', ()), ('INT4H', 'S4', ()), ('TMPDIST4', '<i8', ()),
             ('BEH8', 'S4', ()), ('INT8', 'S4', ()), ('BEH8H', 'S4', ()), ('INT8H', 'S4', ()), ('TMPDIST8', '<i8', ()),
             ('BEH6', 'S4', ()), ('INT6', 'S4', ()), ('BEH6H', 'S4', ()), ('INT6H', 'S4', ()), ('TMPDIST6', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSHLN2(object):
    name = 'PSHLN2'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('DIRECT', '<i8', ()), ('T', '<f8', ()), ('IANAL', '<i8', ()),
             ('BEH3', 'S4', ()), ('INT3', 'S4', ()), ('BEH3H', 'S4', ()), ('INT3H', 'S4', ()), ('BEH4', 'S4', ()),
             ('INT4', 'S4', ()), ('BEH4H', 'S4', ()), ('INT4H', 'S4', ()), ('BEH5', 'S4', ()), ('INT5', 'S4', ()),
             ('BEH5H', 'S4', ()), ('INT5H', 'S4', ()), ('BEH6', 'S4', ()), ('INT6', 'S4', ()), ('BEH6H', 'S4', ()),
             ('INT6H', 'S4', ()), ('BEH8', 'S4', ()), ('INT8', 'S4', ()), ('BEH8H', 'S4', ()), ('INT8H', 'S4', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSLDN1(object):
    name = 'PSLDN1'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('DIRECT', '<i8', ()), ('IANAL', '<i8', ()), ('BEH4', 'S4', ()),
             ('INT4', 'S4', ()), ('BEH4H', 'S4', ()), ('INT4H', 'S4', ()), ('BEH8', 'S4', ()), ('INT8', 'S4', ()),
             ('BEH8H', 'S4', ()), ('INT8H', 'S4', ()), ('BEH10', 'S4', ()), ('INT10', 'S4', ()), ('BEH10H', 'S4', ()),
             ('INT10H', 'S4', ()), ('BEH20', 'S4', ()), ('INT20', 'S4', ()), ('BEH20H', 'S4', ()), ('INT20H', 'S4', ()),
             ('BEH6', 'S4', ()), ('INT6', 'S4', ()), ('BEH6H', 'S4', ()), ('INT6H', 'S4', ()), ('BEH15', 'S4', ()),
             ('INT15', 'S4', ()), ('BEH15H', 'S4', ()), ('INT15H', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSOLID(object):
    name = 'PSOLID'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('CORDM', '<i8', ()), ('IN', '<i8', ()), ('STRESS', '<i8', ()),
             ('ISOP', '<i8', ()), ('FCTN', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PTUBE(object):
    name = 'PTUBE'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('OD', '<f8', ()), ('T', '<f8', ()), ('NSM', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PVISC(object):
    name = 'PVISC'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('CE', '<f8', ()), ('CR', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PWELD(object):
    name = 'PWELD'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('PID', '<i8', ()), ('MID', '<i8', ()), ('D', '<f8', ()), ('CONNBEH', '<i8', ()), ('CONNTYPE', '<i8', ()),
             ('EXTCON', '<i8', ()), ('CONDTYPE', '<i8', ()), ('WELDTYPE', '<i8', ()), ('MINLEN', '<f8', ()),
             ('MAXLEN', '<f8', ()), ('GMCHK', '<i8', ()), ('SPCGS', '<i8', ()), ('CMASS', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SNORM(object):
    name = 'SNORM'
    path = '/NASTRAN/INPUT/PROPERTY'
    dtype = [('GID', '<i8', ()), ('CID', '<i8', ()), ('N1', '<f8', ()), ('N2', '<f8', ()), ('N3', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID(object):
    name = 'GRID'
    path = '/NASTRAN/INPUT/PROPERTY/VCCT'
    dtype = [('GI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SETID(object):
    name = 'SETID'
    path = '/NASTRAN/INPUT/PROPERTY/VCCT'
    dtype = [('SET3ID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THBY(object):
    name = 'THBY'
    path = '/NASTRAN/INPUT/PROPERTY/VCCT'
    dtype = [('G1', '<i8', ()), ('G2', '<i8', ()), ('GINC', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/VCCT'
    dtype = [('ID', '<i8', ()), ('IDCR', '<i8', ()), ('ITYPE', '<i8', ()), ('IGROW', '<i8', ()), ('INCM', '<i8', ()),
             ('METHOD', '<i8', ()), ('TIME', '<f8', ()), ('IACT', '<i8', ()), ('CGI', '<f8', ()), ('GC', '<f8', ()),
             ('GTH', '<f8', ()), ('C', '<f8', ()), ('M', '<f8', ()), ('GMIN', '<f8', ()), ('GC2', '<f8', ()),
             ('GC3', '<f8', ()), ('TABCGI', '<i8', ()), ('TABGC', '<i8', ()), ('TABGTH', '<i8', ()),
             ('TABC', '<i8', ()), ('TABM', '<i8', ()), ('TABGMIN', '<i8', ()), ('TABGC2', '<i8', ()),
             ('TABGC3', '<i8', ()), ('GRID_POS', '<i8', ()), ('GRID_LEN', '<i8', ()), ('THBY_POS', '<i8', ()),
             ('THBY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/VCCT/GRID', '/NASTRAN/INPUT/PROPERTY/VCCT/SETID',
                 '/NASTRAN/INPUT/PROPERTY/VCCT/THBY']


@register_table
class NEVADA(object):
    name = 'NEVADA'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('RENOREF', '<i8', ()), ('RESTART1', '<i8', ()), ('RENORC', '<i8', ()), ('VEGASRC', '<i8', ()),
             ('ECUTOFF', '<i8', ()), ('CONF1', '<f8', ()), ('GRIDCL', '<f8', ()), ('GRIDIT', '<i8', ()),
             ('TIMESC1', '<f8', ()), ('RADKCO1', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SRQ(object):
    name = 'SRQ'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('FLUXSOL', 'S8', ()), ('CONVT', '<f8', ()), ('MAXITER2', '<i8', ()), ('FIJSMM', 'S8', ()),
             ('FILCOFF', '<f8', ()), ('FIJSMT', '<f8', ()), ('FIJSMI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class SRR(object):
    name = 'SRR'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('GEBHART', 'S8', ()), ('CONVTOL', '<f8', ()), ('MAXITER', '<i8', ()), ('FIJSMTM', 'S8', ()),
             ('FIJFILC', '<f8', ()), ('FIJSMTT', '<f8', ()), ('FIJSMTI', '<i8', ()), ('BIJSMTM', 'S8', ()),
             ('BIJFILC', '<f8', ()), ('BIJSMTT', '<f8', ()), ('BIJMAXI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class THERMICA(object):
    name = 'THERMICA'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('SFLUX', '<f8', ()), ('PALBDO', '<f8', ()), ('PBB', '<f8', ()), ('RSTRT', '<i8', ()),
             ('VFART', '<i8', ()), ('RADRAYC', '<i8', ()), ('ORBFLUX', '<i8', ()), ('CONF', '<f8', ()),
             ('TIMESC2', '<f8', ()), ('RADKC2', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TRASYS(object):
    name = 'TRASYS'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('AXIRADM', '<i8', ()), ('AXIAXIM', '<i8', ()), ('AXIANGM', '<i8', ()), ('TIMESC3', '<f8', ()),
             ('RADKC3', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class TSS(object):
    name = 'TSS'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('JUNK', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/PROPERTY/VIEWEX'
    dtype = [('ICAVITY', '<i8', ()), ('RUNINT', '<i8', ()), ('RADKDIS', 'S8', ()), ('ORBITAL', '<i8', ()),
             ('REUSE', '<i8', ()), ('NEVADA', '<i8', ()), ('TSS', '<i8', ()), ('THERMICA', '<i8', ()),
             ('TRASYS', '<i8', ()), ('SRR', '<i8', ()), ('SRQ', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/PROPERTY/VIEWEX/NEVADA', '/NASTRAN/INPUT/PROPERTY/VIEWEX/SRQ',
                 '/NASTRAN/INPUT/PROPERTY/VIEWEX/SRR', '/NASTRAN/INPUT/PROPERTY/VIEWEX/THERMICA',
                 '/NASTRAN/INPUT/PROPERTY/VIEWEX/TRASYS', '/NASTRAN/INPUT/PROPERTY/VIEWEX/TSS']


@register_table
class MKAERO1(object):
    name = 'MKAERO1'
    path = '/NASTRAN/INPUT/TABLE'
    dtype = [('M1', '<f8', ()), ('K1', '<f8', ()), ('U2', '<i8', ()), ('M2', '<f8', ()), ('K2', '<f8', ()),
             ('U3', '<i8', ()), ('M3', '<f8', ()), ('K3', '<f8', ()), ('U4', '<i8', ()), ('M4', '<f8', ()),
             ('K4', '<f8', ()), ('U5', '<i8', ()), ('M5', '<f8', ()), ('K5', '<f8', ()), ('U6', '<i8', ()),
             ('M6', '<f8', ()), ('K6', '<f8', ()), ('U7', '<i8', ()), ('M7', '<f8', ()), ('K7', '<f8', ()),
             ('U8', '<i8', ()), ('M8', '<f8', ()), ('K8', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/INPUT/TABLE/MKAERO2'
    dtype = [('MKS_POS', '<i8', ()), ('MKS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MKS(object):
    name = 'MKS'
    path = '/NASTRAN/INPUT/TABLE/MKAERO2'
    dtype = [('M', '<f8', ()), ('K', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DAMP(object):
    name = 'DAMP'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('F', '<f8', ()), ('G', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABDMP1'
    dtype = [('ID', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABDMP1/DAMP']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D0'
    dtype = [('XI', '<f8', ()), ('YI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D0'
    dtype = [('ITID', '<i8', ()), ('KIND', '<i8', ()), ('EXTRP', '<i8', ()), ('ITIDS', '<i8', ()), ('ITIDB', '<i8', ()),
             ('SM', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABL3D0/XY']


@register_table
class W1(object):
    name = 'W1'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('X1I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W2(object):
    name = 'W2'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('X2I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W3(object):
    name = 'W3'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('X3I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W4(object):
    name = 'W4'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('X4I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class WY(object):
    name = 'WY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('YIJ', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D1'
    dtype = [('ITID', '<i8', ()), ('NV', '<i8', ()), ('KIND1', '<i8', ()), ('KIND2', '<i8', ()), ('KIND3', '<i8', ()),
             ('KIND4', '<i8', ()), ('NW1', '<i8', ()), ('NW2', '<i8', ()), ('NW3', '<i8', ()), ('NW4', '<i8', ()),
             ('EXTRP1', '<i8', ()), ('EXTRP2', '<i8', ()), ('EXTRP3', '<i8', ()), ('EXTRP4', '<i8', ()),
             ('ITIDS1', '<i8', ()), ('ITIDB1', '<i8', ()), ('ITIDS2', '<i8', ()), ('ITIDB2', '<i8', ()),
             ('ITIDS3', '<i8', ()), ('ITIDB3', '<i8', ()), ('ITIDS4', '<i8', ()), ('ITIDB4', '<i8', ()),
             ('SM1', '<i8', ()), ('SM2', '<i8', ()), ('SM3', '<i8', ()), ('SM4', '<i8', ()), ('NWY', '<i8', ()),
             ('W1_POS', '<i8', ()), ('W1_LEN', '<i8', ()), ('W2_POS', '<i8', ()), ('W2_LEN', '<i8', ()),
             ('W3_POS', '<i8', ()), ('W3_LEN', '<i8', ()), ('W4_POS', '<i8', ()), ('W4_LEN', '<i8', ()),
             ('WY_POS', '<i8', ()), ('WY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABL3D1/W1', '/NASTRAN/INPUT/TABLE/TABL3D1/W2',
                 '/NASTRAN/INPUT/TABLE/TABL3D1/W3', '/NASTRAN/INPUT/TABLE/TABL3D1/W4',
                 '/NASTRAN/INPUT/TABLE/TABL3D1/WY']


@register_table
class W1(object):
    name = 'W1'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('X1I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W2(object):
    name = 'W2'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('X2I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W3(object):
    name = 'W3'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('X3I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class W4(object):
    name = 'W4'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('X4I', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class WY(object):
    name = 'WY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('YIJ', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D2'
    dtype = [('ITID', '<i8', ()), ('NV', '<i8', ()), ('KIND1', '<i8', ()), ('KIND2', '<i8', ()), ('KIND3', '<i8', ()),
             ('KIND4', '<i8', ()), ('NW1', '<i8', ()), ('NW2', '<i8', ()), ('NW3', '<i8', ()), ('NW4', '<i8', ()),
             ('EXTRP1', '<i8', ()), ('EXTRP2', '<i8', ()), ('EXTRP3', '<i8', ()), ('EXTRP4', '<i8', ()),
             ('ITIDS1', '<i8', ()), ('ITIDB1', '<i8', ()), ('ITIDS2', '<i8', ()), ('ITIDB2', '<i8', ()),
             ('ITIDS3', '<i8', ()), ('ITIDB3', '<i8', ()), ('ITIDS4', '<i8', ()), ('ITIDB4', '<i8', ()),
             ('SM1', '<i8', ()), ('SM2', '<i8', ()), ('SM3', '<i8', ()), ('SM4', '<i8', ()), ('NWY', '<i8', ()),
             ('W1_POS', '<i8', ()), ('W1_LEN', '<i8', ()), ('W2_POS', '<i8', ()), ('W2_LEN', '<i8', ()),
             ('W3_POS', '<i8', ()), ('W3_LEN', '<i8', ()), ('W4_POS', '<i8', ()), ('W4_LEN', '<i8', ()),
             ('WY_POS', '<i8', ()), ('WY_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABL3D2/W1', '/NASTRAN/INPUT/TABLE/TABL3D2/W2',
                 '/NASTRAN/INPUT/TABLE/TABL3D2/W3', '/NASTRAN/INPUT/TABLE/TABL3D2/W4',
                 '/NASTRAN/INPUT/TABLE/TABL3D2/WY']


@register_table
class F(object):
    name = 'F'
    path = '/NASTRAN/INPUT/TABLE/TABL3D3'
    dtype = [('FORMULA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABL3D3'
    dtype = [('ITID', '<i8', ()), ('NV', '<i8', ()), ('KIND1', '<i8', ()), ('KIND2', '<i8', ()), ('KIND3', '<i8', ()),
             ('KIND4', '<i8', ()), ('NW1', '<i8', ()), ('NW2', '<i8', ()), ('NW3', '<i8', ()), ('NW4', '<i8', ()),
             ('EXTRP1', '<i8', ()), ('EXTRP2', '<i8', ()), ('EXTRP3', '<i8', ()), ('EXTRP4', '<i8', ()),
             ('ITIDS1', '<i8', ()), ('ITIDB1', '<i8', ()), ('ITIDS2', '<i8', ()), ('ITIDB2', '<i8', ()),
             ('ITIDS3', '<i8', ()), ('ITIDB3', '<i8', ()), ('ITIDS4', '<i8', ()), ('ITIDB4', '<i8', ()),
             ('SM1', '<i8', ()), ('SM2', '<i8', ()), ('SM3', '<i8', ()), ('SM4', '<i8', ()), ('POS', '<i8', ()),
             ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABL3D3/F']


@register_table
class ENTRIES(object):
    name = 'ENTRIES'
    path = '/NASTRAN/INPUT/TABLE/TABLE3D'
    dtype = [('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()), ('FI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLE3D'
    dtype = [('ID', '<i8', ()), ('X0', '<f8', ()), ('Y0', '<f8', ()), ('Z0', '<f8', ()), ('F0', '<f8', ()),
             ('ENTRIES_POS', '<i8', ()), ('ENTRIES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLE3D/ENTRIES']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLED1'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLED1'
    dtype = [('ID', '<i8', ()), ('CODEX', '<i8', ()), ('CODEY', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLED1/XY']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLED2'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLED2'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLED2/XY']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLED3'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLED3'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLED3/XY']


@register_table
class COEF(object):
    name = 'COEF'
    path = '/NASTRAN/INPUT/TABLE/TABLED4'
    dtype = [('A', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLED4'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('X4', '<f8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLED4/COEF']


@register_table
class ENTRIES(object):
    name = 'ENTRIES'
    path = '/NASTRAN/INPUT/TABLE/TABLED5'
    dtype = [('XI', '<f8', ()), ('TID', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLED5'
    dtype = [('TID', '<i8', ()), ('ENTRIES_POS', '<i8', ()), ('ENTRIES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLED5/ENTRIES']


@register_table
class ENTRIES(object):
    name = 'ENTRIES'
    path = '/NASTRAN/INPUT/TABLE/TABLEH1'
    dtype = [('YI', '<f8', ()), ('FI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEH1'
    dtype = [('TIB', '<i8', ()), ('ENTRIES_POS', '<i8', ()), ('ENTRIES_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEH1/ENTRIES']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM1'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM1'
    dtype = [('ID', '<i8', ()), ('CODEX', '<i8', ()), ('CODEY', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEM1/XY']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM2'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM2'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEM2/XY']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM3'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM3'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEM3/XY']


@register_table
class COEF(object):
    name = 'COEF'
    path = '/NASTRAN/INPUT/TABLE/TABLEM4'
    dtype = [('A', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEM4'
    dtype = [('ID', '<i8', ()), ('X1', '<f8', ()), ('X2', '<f8', ()), ('X3', '<f8', ()), ('X4', '<f8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEM4/COEF']


@register_table
class XY(object):
    name = 'XY'
    path = '/NASTRAN/INPUT/TABLE/TABLES1'
    dtype = [('X', '<f8', ()), ('Y', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLES1'
    dtype = [('ID', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLES1/XY']


@register_table
class T(object):
    name = 'T'
    path = '/NASTRAN/INPUT/TABLE/TABLEST'
    dtype = [('TI', '<f8', ()), ('TIDI', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLEST'
    dtype = [('ID', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLEST/T']


@register_table
class YIS(object):
    name = 'YIS'
    path = '/NASTRAN/INPUT/TABLE/TABLFTG'
    dtype = [('YI', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABLFTG'
    dtype = [('TID', '<i8', ()), ('YIS_POS', '<i8', ()), ('YIS_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABLFTG/YIS']


@register_table
class F(object):
    name = 'F'
    path = '/NASTRAN/INPUT/TABLE/TABRND1'
    dtype = [('F', '<f8', ()), ('G', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABRND1'
    dtype = [('ID', '<i8', ()), ('CODEX', '<i8', ()), ('CODEY', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABRND1/F']


@register_table
class TABRNDG(object):
    name = 'TABRNDG'
    path = '/NASTRAN/INPUT/TABLE'
    dtype = [('TID', '<i8', ()), ('TYPE', '<i8', ()), ('LU', '<f8', ()), ('WG', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CRITERIA(object):
    name = 'CRITERIA'
    path = '/NASTRAN/INPUT/TABLE/TABSCTL'
    dtype = [('ICRIT', 'S8', ()), ('SET3', '<i8', ()), ('YT', '<f8', (4,)), ('MT', '<f8', (4,))]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/TABLE/TABSCTL'
    dtype = [('ID', '<i8', ()), ('CRITERIA_POS', '<i8', ()), ('CRITERIA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/TABLE/TABSCTL/CRITERIA']


@register_table
class CHARS(object):
    name = 'CHARS'
    path = '/NASTRAN/INPUT/UDS/ELEMUDS'
    dtype = [('CDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class INTS(object):
    name = 'INTS'
    path = '/NASTRAN/INPUT/UDS/ELEMUDS'
    dtype = [('IDATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class REALS(object):
    name = 'REALS'
    path = '/NASTRAN/INPUT/UDS/ELEMUDS'
    dtype = [('RDATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/UDS/ELEMUDS'
    dtype = [('NWD', '<i8', ()), ('PID', '<i8', ()), ('PTYPE', 'S8', ()), ('GROUP', 'S8', ()), ('UNAME', 'S8', ()),
             ('DEPEN', 'S8', ()), ('NINTG', '<i8', ()), ('NREAL', '<i8', ()), ('NCHAR', '<i8', ()),
             ('INTS_POS', '<i8', ()), ('REALS_POS', '<i8', ()), ('CHARS_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/UDS/ELEMUDS/CHARS', '/NASTRAN/INPUT/UDS/ELEMUDS/INTS',
                 '/NASTRAN/INPUT/UDS/ELEMUDS/REALS']


@register_table
class CHARS(object):
    name = 'CHARS'
    path = '/NASTRAN/INPUT/UDS/ENTUDS'
    dtype = [('CDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class INTS(object):
    name = 'INTS'
    path = '/NASTRAN/INPUT/UDS/ENTUDS'
    dtype = [('IDATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class REALS(object):
    name = 'REALS'
    path = '/NASTRAN/INPUT/UDS/ENTUDS'
    dtype = [('RDATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/UDS/ENTUDS'
    dtype = [('NWD', '<i8', ()), ('ENTID', '<i8', ()), ('ENTPNT', 'S8', ()), ('ENTGN', 'S8', ()), ('NINTG', '<i8', ()),
             ('NREAL', '<i8', ()), ('NCHAR', '<i8', ()), ('INTS_POS', '<i8', ()), ('REALS_POS', '<i8', ()),
             ('CHARS_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/UDS/ENTUDS/CHARS', '/NASTRAN/INPUT/UDS/ENTUDS/INTS', '/NASTRAN/INPUT/UDS/ENTUDS/REALS']


@register_table
class CHARS(object):
    name = 'CHARS'
    path = '/NASTRAN/INPUT/UDS/GENUDS'
    dtype = [('CDATA', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class INTS(object):
    name = 'INTS'
    path = '/NASTRAN/INPUT/UDS/GENUDS'
    dtype = [('IDATA', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class REALS(object):
    name = 'REALS'
    path = '/NASTRAN/INPUT/UDS/GENUDS'
    dtype = [('RDATA', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/UDS/GENUDS'
    dtype = [('NWD', '<i8', ()), ('SRV', 'S8', ()), ('NINTG', '<i8', ()), ('NREAL', '<i8', ()), ('NCHAR', '<i8', ()),
             ('INTS_POS', '<i8', ()), ('REALS_POS', '<i8', ()), ('CHARS_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/UDS/GENUDS/CHARS', '/NASTRAN/INPUT/UDS/GENUDS/INTS', '/NASTRAN/INPUT/UDS/GENUDS/REALS']


@register_table
class CASES(object):
    name = 'CASES'
    path = '/NASTRAN/INPUT/UDS/THPAD'
    dtype = [('RPM', '<f8', ()), ('FX', '<f8', ()), ('FY', '<f8', ()), ('PGNU', '<f8', ()), ('RELAX', '<f8', ()),
             ('OFLOW', '<f8', ()), ('DEREL', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class PADS(object):
    name = 'PADS'
    path = '/NASTRAN/INPUT/UDS/THPAD'
    dtype = [('ARC', '<f8', ()), ('OFFSET', '<f8', ()), ('PRELOAD', '<f8', ()), ('PVANG', '<f8', ()), ('IP', '<f8', ()),
             ('KP', '<f8', ()), ('MP', '<f8', ()), ('IT', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/UDS/THPAD'
    dtype = [('NWD', '<i8', ()), ('RID', '<i8', ()), ('TLIMIT', '<f8', ()), ('TITLE1', 'S64', ()),
             ('TITLE2', 'S64', ()), ('RADIUS', '<f8', ()), ('CLEAR', '<f8', ()), ('OUTR', '<f8', ()),
             ('LENGTH', '<f8', ()), ('E', '<f8', ()), ('ALPHP', '<f8', ()), ('ALPHJ', '<f8', ()), ('ALPHS', '<f8', ()),
             ('NPADS', '<i8', ()), ('NEL', '<i8', ()), ('IECC', '<i8', ()), ('IPLOT', '<i8', ()), ('ISUMRY', '<i8', ()),
             ('PRNFULL', '<i8', ()), ('PADS_POS', '<i8', ()), ('KFILM', '<f8', ()), ('KPAD', '<f8', ()),
             ('TBACK', '<f8', ()), ('TJF', '<f8', ()), ('PSUMP', '<f8', ()), ('TIN', '<f8', ()), ('KCAV', '<f8', ()),
             ('CCAV', '<f8', ()), ('DENSITY', '<f8', ()), ('SPEC', '<f8', ()), ('TA', '<f8', ()), ('MUA', '<f8', ()),
             ('TB', '<f8', ()), ('MUB', '<f8', ()), ('ESUMP', '<f8', ()), ('TMANU', '<f8', ()), ('ERROR', '<f8', ()),
             ('XG', '<f8', ()), ('YG', '<f8', ()), ('NAX', '<f8', ()), ('FACTOR', '<f8', ()), ('XFACT', '<f8', ()),
             ('ITJ', '<i8', ()), ('IBC', '<i8', ()), ('ITUR', '<i8', ()), ('IDIM', '<i8', ()), ('ICOND', '<i8', ()),
             ('ITB', '<i8', ()), ('ITG', '<i8', ()), ('IOUT1', '<i8', ()), ('IOUT2', '<i8', ()), ('ICAV', '<i8', ()),
             ('ICROSS', '<i8', ()), ('IDEF', '<i8', ()), ('IFLEX', '<i8', ()), ('IUN', '<i8', ()), ('ISUMP', '<i8', ()),
             ('KTHETA', '<f8', ()), ('DEL0', '<f8', ()), ('HOTOVER', '<f8', ()), ('COLDOVER', '<f8', ()),
             ('NCASE', '<i8', ()), ('CASES_POS', '<i8', ()), ('MAXC', '<i8', ()), ('IALPH', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/UDS/THPAD/CASES', '/NASTRAN/INPUT/UDS/THPAD/PADS']


@register_table
class SVNAMES(object):
    name = 'SVNAMES'
    path = '/NASTRAN/INPUT/UDS/UDSESV'
    dtype = [('SVI', 'S8', ()), ('SVI_NAME', 'S8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/INPUT/UDS/UDSESV'
    dtype = [('NSTATS', '<i8', ()), ('NAME_CNT', '<i8', ()), ('SVNAMES_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/INPUT/UDS/UDSESV/SVNAMES']


@register_table
class AFPM(object):
    name = 'AFPM'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('PR', '<f8', ()), ('PI', '<f8', ()), ('IN', '<f8', ()), ('IX', '<f8', ()),
             ('IY', '<f8', ()), ('IZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ERP(object):
    name = 'ERP'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('PANEL_NAME', 'S8', ()), ('ERPVAL', '<f8', ()), ('ERPFRA', '<f8', ()), ('ERPdB', '<f8', ()),
             ('AREA', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/ACOUSTIC/ERPE'
    dtype = [('IDENT', '<i8', ()), ('ERPMAX', '<f8', ()), ('FLTRLOC', '<i8', ()), ('FILTER', '<f8', ()),
             ('PANEL', 'S8', ()), ('ERP', '<f8', ()), ('AREA', '<f8', ()), ('SORTKEY', '<i8', ()),
             ('ASORDS', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class POWER(object):
    name = 'POWER'
    path = '/NASTRAN/RESULT/ACOUSTIC/ERPE'
    dtype = [('EID', '<i8', ()), ('ERPVAL', '<f8', ()), ('ERPFRA', '<f8', ()), ('ERPdB', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/ACOUSTIC/ERPM'
    dtype = [('IDENT', '<i8', ()), ('ERPMAX', '<f8', ()), ('FLTRLOC', '<i8', ()), ('FILTER', '<f8', ()),
             ('PANEL', 'S8', ()), ('ERP', '<f8', ()), ('AREA', '<f8', ()), ('SORTKEY', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class POWER(object):
    name = 'POWER'
    path = '/NASTRAN/RESULT/ACOUSTIC/ERPM'
    dtype = [('EIGFREQ', '<f8', ()), ('ERPVAL', '<f8', ()), ('ERPFRA', '<f8', ()), ('ERPdB', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXPR(object):
    name = 'HEXPR'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('EID', '<i8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()), ('AZ', '<f8', ()), ('VX', '<f8', ()),
             ('VY', '<f8', ()), ('VZ', '<f8', ()), ('PRESSURE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXPRINT(object):
    name = 'HEXPRINT'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('EID', '<i8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()), ('AZ', '<f8', ()), ('VX', '<f8', ()),
             ('VY', '<f8', ()), ('VZ', '<f8', ()), ('PRESSURE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXPRINT_CPLX(object):
    name = 'HEXPRINT_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('EID', '<i8', ()), ('AXR', '<f8', ()), ('AYR', '<f8', ()), ('AZR', '<f8', ()), ('VXR', '<f8', ()),
             ('VYR', '<f8', ()), ('VZR', '<f8', ()), ('PRESSURE', '<f8', ()), ('AXI', '<f8', ()), ('AYI', '<f8', ()),
             ('AZI', '<f8', ()), ('VXI', '<f8', ()), ('VYI', '<f8', ()), ('VZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXPR_CPLX(object):
    name = 'HEXPR_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('EID', '<i8', ()), ('AXR', '<f8', ()), ('AYR', '<f8', ()), ('AZR', '<f8', ()), ('VXR', '<f8', ()),
             ('VYR', '<f8', ()), ('VZR', '<f8', ()), ('PRESSURE', '<f8', ()), ('AXI', '<f8', ()), ('AYI', '<f8', ()),
             ('AZI', '<f8', ()), ('VXI', '<f8', ()), ('VYI', '<f8', ()), ('VZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class INTENSITY(object):
    name = 'INTENSITY'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('IN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLUID_MODAL(object):
    name = 'FLUID_MODAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF/STRUCTURE_MODAL'
    subtables = []


@register_table
class GRID(object):
    name = 'GRID'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = [('GRID', '<i8', ()), ('TYPE', 'S4', ()), ('DT1R', '<f8', ()), ('DT2R', '<f8', ()), ('DT3R', '<f8', ()),
             ('DR1R', '<f8', ()), ('DR2R', '<f8', ()), ('DR3R', '<f8', ()), ('DT1I', '<f8', ()), ('DT2I', '<f8', ()),
             ('DT3I', '<f8', ()), ('DR1I', '<f8', ()), ('DR2I', '<f8', ()), ('DR3I', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = [('IDENT', '<i8', ()), ('GRID', '<i8', ()), ('C', '<i8', ()), ('FILTER', '<f8', ()), ('RR', '<f8', ()),
             ('RI', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PANEL(object):
    name = 'PANEL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = [('PID', '<i8', ()), ('PANEL', 'S8', ()), ('PR', '<f8', ()), ('PI', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PANEL_MODAL(object):
    name = 'PANEL_MODAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF/STRUCTURE_MODAL'
    subtables = []


@register_table
class STRUCTURE_MODAL(object):
    name = 'STRUCTURE_MODAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = [('MODE', '<i8', ()), ('RFREQ', '<f8', ()), ('PR', '<f8', ()), ('PI', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class STRUCTURE_PANEL(object):
    name = 'STRUCTURE_PANEL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF/PANEL'
    subtables = []


@register_table
class STRUCTURE_REAL(object):
    name = 'STRUCTURE_REAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF'
    dtype = [('MODE', '<i8', ()), ('RFREQ', '<f8', ()), ('PR', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLUID_FORCING(object):
    name = 'FLUID_FORCING'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class FLUID_NATURAL(object):
    name = 'FLUID_NATURAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = [('FREQUENCY', '<f8', ()), ('TYPE', '<i8', ()), ('MPFR', '<f8', ()), ('MPFI', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_PANEL_FORCING(object):
    name = 'GRID_PANEL_FORCING'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class GRID_PANEL_NATURAL(object):
    name = 'GRID_PANEL_NATURAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = [('IDENT', '<i8', ()), ('ID', '<i8', ()), ('FREQ', '<f8', ()), ('FRQI', '<i8', ()), ('PANEL', 'S8', ()),
             ('FASETDOF', '<i8', ()), ('SUBCASE', '<i8', ()), ('PGRIDID', '<i8', ()), ('SPGID', '<i8', ()),
             ('PANELI', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LOAD_FORCING(object):
    name = 'LOAD_FORCING'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class LOAD_POINT(object):
    name = 'LOAD_POINT'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = [('ID', '<i8', ()), ('TYPE', '<i8', ()), ('MPFR', '<i8', ()), ('MPFI', '<i8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PANEL_FORCING(object):
    name = 'PANEL_FORCING'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class PANEL_NATURAL(object):
    name = 'PANEL_NATURAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class STRUCTURE_FORCING(object):
    name = 'STRUCTURE_FORCING'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class STRUCTURE_NATURAL(object):
    name = 'STRUCTURE_NATURAL'
    path = '/NASTRAN/RESULT/ACOUSTIC/MPF2'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/MPF2/FLUID_NATURAL'
    subtables = []


@register_table
class PENPR(object):
    name = 'PENPR'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPR'
    subtables = []


@register_table
class PENPRINT(object):
    name = 'PENPRINT'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPRINT'
    subtables = []


@register_table
class PENPRINT_CPLX(object):
    name = 'PENPRINT_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPRINT_CPLX'
    subtables = []


@register_table
class PENPR_CPLX(object):
    name = 'PENPR_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPR_CPLX'
    subtables = []


@register_table
class POWER(object):
    name = 'POWER'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('PANEL_ID', '<i8', ()), ('PANEL_NAME', 'S8', ()), ('POWER', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class POWER_FPM(object):
    name = 'POWER_FPM'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('POWER', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRESSURE(object):
    name = 'PRESSURE'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('P', '<f8', ()), ('PRMS', '<f8', ()), ('DB', '<f8', ()), ('DBA', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRESSURE_CPLX(object):
    name = 'PRESSURE_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('PR', '<f8', ()), ('PRMSR', '<f8', ()), ('DBR', '<f8', ()), ('DBAR', '<f8', ()),
             ('PI', '<f8', ()), ('PRMSI', '<f8', ()), ('DBI', '<f8', ()), ('DBAI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETPR(object):
    name = 'TETPR'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPR'
    subtables = []


@register_table
class TETPRINT(object):
    name = 'TETPRINT'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPRINT'
    subtables = []


@register_table
class TETPRINT_CPLX(object):
    name = 'TETPRINT_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPRINT_CPLX'
    subtables = []


@register_table
class TETPR_CPLX(object):
    name = 'TETPR_CPLX'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ACOUSTIC/HEXPR_CPLX'
    subtables = []


@register_table
class VELOCITY_FPM(object):
    name = 'VELOCITY_FPM'
    path = '/NASTRAN/RESULT/ACOUSTIC'
    dtype = [('ID', '<i8', ()), ('VT1R', '<f8', ()), ('VT2R', '<f8', ()), ('VT3R', '<f8', ()), ('VT1I', '<f8', ()),
             ('VT2I', '<f8', ()), ('VT3I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLEXIBLE(object):
    name = 'FLEXIBLE'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('FBCID', '<i8', ()), ('SBCID', '<i8', ()), ('TBCID', '<i8', ()),
             ('CONSTAT', '<i8', ()), ('FCONN', '<f8', ()), ('FCONF', '<f8', ()), ('SNORM', '<f8', ()),
             ('SFRIC1', '<f8', ()), ('SFRIC2', '<f8', ()), ('FNX', '<f8', ()), ('FNY', '<f8', ()), ('FNZ', '<f8', ()),
             ('FFX', '<f8', ()), ('FFY', '<f8', ()), ('FFZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GLUED_NORM(object):
    name = 'GLUED_NORM'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GLUED_NORM_CPLX(object):
    name = 'GLUED_NORM_CPLX'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GLUED_TANG(object):
    name = 'GLUED_TANG'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GLUED_TANG_CPLX(object):
    name = 'GLUED_TANG_CPLX'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RIGID(object):
    name = 'RIGID'
    path = '/NASTRAN/RESULT/CONTACT'
    dtype = [('ID', '<i8', ()), ('XMAT1', '<f8', ()), ('XMAT2', '<f8', ()), ('XMAT3', '<f8', ()), ('XMAT4', '<f8', ()),
             ('XMAT5', '<f8', ()), ('XMAT6', '<f8', ()), ('XMAT7', '<f8', ()), ('XMAT8', '<f8', ()),
             ('XMAT9', '<f8', ()), ('XMAT10', '<f8', ()), ('XMAT11', '<f8', ()), ('XMAT12', '<f8', ()),
             ('XMAT13', '<f8', ()), ('XMAT14', '<f8', ()), ('XMAT15', '<f8', ()), ('XMAT16', '<f8', ()),
             ('RFX', '<f8', ()), ('RFY', '<f8', ()), ('RFZ', '<f8', ()), ('RMX', '<f8', ()), ('RMY', '<f8', ()),
             ('RMZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DOMAINS(object):
    name = 'DOMAINS'
    path = '/NASTRAN/RESULT'
    dtype = [('ID', '<i8', ()), ('SUBCASE', '<i8', ()), ('STEP', '<i8', ()), ('ANALYSIS', '<i8', ()),
             ('TIME_FREQ_EIGR', '<f8', ()), ('EIGI', '<f8', ()), ('MODE', '<i8', ()), ('DESIGN_CYCLE', '<i8', ()),
             ('RANDOM', '<i8', ()), ('SE', '<i8', ()), ('AFPM', '<i8', ()), ('TRMC', '<i8', ()),
             ('INSTANCE', '<i8', ()), ('MODULE', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR(object):
    name = 'BAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('BM1A', '<f8', ()), ('BM2A', '<f8', ()), ('BM1B', '<f8', ()), ('BM2B', '<f8', ()),
             ('TS1', '<f8', ()), ('TS2', '<f8', ()), ('AF', '<f8', ()), ('TRQ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BARS(object):
    name = 'BARS'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('SD', '<f8', ()), ('BM1', '<f8', ()), ('BM2', '<f8', ()), ('TS1', '<f8', ()),
             ('TS2', '<f8', ()), ('AF', '<f8', ()), ('TRQ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BARS_CPLX(object):
    name = 'BARS_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('SD', '<f8', ()), ('BM1R', '<f8', ()), ('BM2R', '<f8', ()), ('TS1R', '<f8', ()),
             ('TS2R', '<f8', ()), ('AFR', '<f8', ()), ('TRQR', '<f8', ()), ('BM1I', '<f8', ()), ('BM2I', '<f8', ()),
             ('TS1I', '<f8', ()), ('TS2I', '<f8', ()), ('AFI', '<f8', ()), ('TRQI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR_CPLX(object):
    name = 'BAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('BM1AR', '<f8', ()), ('BM2AR', '<f8', ()), ('BM1BR', '<f8', ()), ('BM2BR', '<f8', ()),
             ('TS1R', '<f8', ()), ('TS2R', '<f8', ()), ('AFR', '<f8', ()), ('TRQR', '<f8', ()), ('BM1AI', '<f8', ()),
             ('BM2AI', '<f8', ()), ('BM1BI', '<f8', ()), ('BM2BI', '<f8', ()), ('TS1I', '<f8', ()), ('TS2I', '<f8', ()),
             ('AFI', '<f8', ()), ('TRQI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM(object):
    name = 'BEAM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('SD', '<f8', ()), ('BM1', '<f8', ()), ('BM2', '<f8', ()),
                                   ('TS1', '<f8', ()), ('TS2', '<f8', ()), ('AF', '<f8', ()), ('TTRQ', '<f8', ()),
                                   ('WTRQ', '<f8', ())], (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM3(object):
    name = 'BEAM3'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('MY', '<f8', ()), ('MZ', '<f8', ()), ('QY', '<f8', ()),
                                   ('QZ', '<f8', ()), ('NX', '<f8', ()), ('TX', '<f8', ()), ('VB', '<f8', ()),
                                   ('MB', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM3_CPLX(object):
    name = 'BEAM3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('MYR', '<f8', ()), ('MZR', '<f8', ()), ('QYR', '<f8', ()),
                                   ('QZR', '<f8', ()), ('NXR', '<f8', ()), ('TXR', '<f8', ()), ('VBR', '<f8', ()),
                                   ('MBR', '<f8', ()), ('MYI', '<f8', ()), ('MZI', '<f8', ()), ('QYI', '<f8', ()),
                                   ('QZI', '<f8', ()), ('NXI', '<f8', ()), ('TXI', '<f8', ()), ('VBI', '<f8', ()),
                                   ('MBI', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_CPLX(object):
    name = 'BEAM_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('SD', '<f8', ()), ('BM1R', '<f8', ()), ('BM2R', '<f8', ()),
                                   ('TS1R', '<f8', ()), ('TS2R', '<f8', ()), ('AFR', '<f8', ()), ('TTRQR', '<f8', ()),
                                   ('WTRQR', '<f8', ()), ('BM1I', '<f8', ()), ('BM2I', '<f8', ()), ('TS1I', '<f8', ()),
                                   ('TS2I', '<f8', ()), ('AFI', '<f8', ()), ('TTRQI', '<f8', ()), ('WTRQI', '<f8', ())],
                                  (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEND(object):
    name = 'BEND'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('BM1', '<f8', ()), ('BM2', '<f8', ()), ('TS1', '<f8', ()),
                                   ('TS2', '<f8', ()), ('AF', '<f8', ()), ('TRQ', '<f8', ())], (2,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEND_CPLX(object):
    name = 'BEND_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('GRID', '<i8', ()), ('BM1R', '<f8', ()), ('BM2R', '<f8', ()), ('TS1R', '<f8', ()),
                                   ('TS2R', '<f8', ()), ('AFR', '<f8', ()), ('TRQR', '<f8', ()), ('BM1I', '<f8', ()),
                                   ('BM2I', '<f8', ()), ('TS1I', '<f8', ()), ('TS2I', '<f8', ()), ('AFI', '<f8', ()),
                                   ('TRQI', '<f8', ())], (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH(object):
    name = 'BUSH'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FX', '<f8', ()), ('FY', '<f8', ()), ('FZ', '<f8', ()), ('MX', '<f8', ()),
             ('MY', '<f8', ()), ('MZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH1D(object):
    name = 'BUSH1D'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FE', '<f8', ()), ('UE', '<f8', ()), ('VE', '<f8', ()), ('AS', '<f8', ()),
             ('AE', '<f8', ()), ('EP', '<f8', ()), ('FAIL', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH_CPLX(object):
    name = 'BUSH_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FXR', '<f8', ()), ('FYR', '<f8', ()), ('FZR', '<f8', ()), ('MXR', '<f8', ()),
             ('MYR', '<f8', ()), ('MZR', '<f8', ()), ('FXI', '<f8', ()), ('FYI', '<f8', ()), ('FZI', '<f8', ()),
             ('MXI', '<f8', ()), ('MYI', '<f8', ()), ('MZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONE(object):
    name = 'CONE'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('HOPA', '<f8', ()), ('BMU', '<f8', ()), ('BMV', '<f8', ()), ('TM', '<f8', ()),
             ('SU', '<f8', ()), ('SV', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD(object):
    name = 'CONROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD'
    subtables = []


@register_table
class CONROD_CPLX(object):
    name = 'CONROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD_CPLX'
    subtables = []


@register_table
class CONV(object):
    name = 'CONV'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('HFLOW', '<f8', ()), ('CNTLNODE', '<f8', ()), ('COEF', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONVM(object):
    name = 'CONVM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/CONV'
    subtables = []


@register_table
class DAMP1(object):
    name = 'DAMP1'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class DAMP1_CPLX(object):
    name = 'DAMP1_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class DAMP2(object):
    name = 'DAMP2'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class DAMP2_CPLX(object):
    name = 'DAMP2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class DAMP3(object):
    name = 'DAMP3'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class DAMP3_CPLX(object):
    name = 'DAMP3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class DAMP4(object):
    name = 'DAMP4'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class DAMP4_CPLX(object):
    name = 'DAMP4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS1(object):
    name = 'ELAS1'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('F', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS1_CPLX(object):
    name = 'ELAS1_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FR', '<f8', ()), ('FI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS2(object):
    name = 'ELAS2'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class ELAS2_CPLX(object):
    name = 'ELAS2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS3(object):
    name = 'ELAS3'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class ELAS3_CPLX(object):
    name = 'ELAS3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS4(object):
    name = 'ELAS4'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1'
    subtables = []


@register_table
class ELAS4_CPLX(object):
    name = 'ELAS4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX'
    subtables = []


@register_table
class FAST(object):
    name = 'FAST'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BUSH'
    subtables = []


@register_table
class FAST_CPLX(object):
    name = 'FAST_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BUSH_CPLX'
    subtables = []


@register_table
class GAP(object):
    name = 'GAP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FX', '<f8', ()), ('SFY', '<f8', ()), ('SFZ', '<f8', ()), ('U', '<f8', ()),
             ('V', '<f8', ()), ('W', '<f8', ()), ('SV', '<f8', ()), ('SW', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRAD_FLUX(object):
    name = 'GRAD_FLUX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('XGRAD', '<f8', ()), ('YGRAD', '<f8', ()), ('ZGRAD', '<f8', ()), ('XFLUX', '<f8', ()),
             ('YFLUX', '<f8', ()), ('ZFLUX', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HBDYE(object):
    name = 'HBDYE'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FAPPLIED', '<f8', ()), ('FREECONV', '<f8', ()), ('FORCECON', '<f8', ()),
             ('FRAD', '<f8', ()), ('FTOTAL', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HBDYG(object):
    name = 'HBDYG'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/HBDYE'
    subtables = []


@register_table
class HBDYP(object):
    name = 'HBDYP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/HBDYE'
    subtables = []


@register_table
class QUAD4(object):
    name = 'QUAD4'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('MX', '<f8', ()), ('MY', '<f8', ()), ('MXY', '<f8', ()), ('BMX', '<f8', ()),
                                   ('BMY', '<f8', ()), ('BMXY', '<f8', ()), ('TX', '<f8', ()), ('TY', '<f8', ())], ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_CN(object):
    name = 'QUAD4_CN'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('FORCE',
                                                      [('GRID', '<i8', ()), ('MX', '<f8', ()), ('MY', '<f8', ()),
                                                       ('MXY', '<f8', ()), ('BMX', '<f8', ()), ('BMY', '<f8', ()),
                                                       ('BMXY', '<f8', ()), ('TX', '<f8', ()), ('TY', '<f8', ())],
                                                      (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_CN_CPLX(object):
    name = 'QUAD4_CN_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('FORCE',
                                                      [('GRID', '<i8', ()), ('MXR', '<f8', ()), ('MYR', '<f8', ()),
                                                       ('MXYR', '<f8', ()), ('BMXR', '<f8', ()), ('BMYR', '<f8', ()),
                                                       ('BMXYR', '<f8', ()), ('TXR', '<f8', ()), ('TYR', '<f8', ()),
                                                       ('MXI', '<f8', ()), ('MYI', '<f8', ()), ('MXYI', '<f8', ()),
                                                       ('BMXI', '<f8', ()), ('BMYI', '<f8', ()), ('BMXYI', '<f8', ()),
                                                       ('TXI', '<f8', ()), ('TYI', '<f8', ())], (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_COMP(object):
    name = 'QUAD4_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('THEORY', 'S8', ()), ('LAMID', '<i8', ()), ('FP', '<f8', ()), ('FM', '<f8', ()),
             ('FB', '<f8', ()), ('FMAX', '<f8', ()), ('FFLAG', 'S4', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_CPLX(object):
    name = 'QUAD4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FORCE',
                                  [('MXR', '<f8', ()), ('MYR', '<f8', ()), ('MXYR', '<f8', ()), ('BMXR', '<f8', ()),
                                   ('BMYR', '<f8', ()), ('BMXYR', '<f8', ()), ('TXR', '<f8', ()), ('TYR', '<f8', ()),
                                   ('MXI', '<f8', ()), ('MYI', '<f8', ()), ('MXYI', '<f8', ()), ('BMXI', '<f8', ()),
                                   ('BMYI', '<f8', ()), ('BMXYI', '<f8', ()), ('TXI', '<f8', ()), ('TYI', '<f8', ())],
                                  ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8(object):
    name = 'QUAD8'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN'
    subtables = []


@register_table
class QUAD8_COMP(object):
    name = 'QUAD8_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_COMP'
    subtables = []


@register_table
class QUAD8_CPLX(object):
    name = 'QUAD8_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN_CPLX'
    subtables = []


@register_table
class QUADR(object):
    name = 'QUADR'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN'
    subtables = []


@register_table
class QUADR_CPLX(object):
    name = 'QUADR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN_CPLX'
    subtables = []


@register_table
class RAC2D(object):
    name = 'RAC2D'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('F', '<f8', (9,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RAC2D_CPLX(object):
    name = 'RAC2D_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('FR', '<f8', (9,)), ('FI', '<f8', (9,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RAC3D(object):
    name = 'RAC3D'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/RAC2D'
    subtables = []


@register_table
class RAC3D_CPLX(object):
    name = 'RAC3D_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/RAC2D_CPLX'
    subtables = []


@register_table
class RADBC(object):
    name = 'RADBC'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('HFLOW', '<f8', ()), ('CNTLNODE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RADINT(object):
    name = 'RADINT'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/RADBC'
    subtables = []


@register_table
class ROD(object):
    name = 'ROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('AF', '<f8', ()), ('TRQ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROD_CPLX(object):
    name = 'ROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('AFR', '<f8', ()), ('AFI', '<f8', ()), ('TRQR', '<f8', ()), ('TRQI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SEAM(object):
    name = 'SEAM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR'
    subtables = []


@register_table
class SEAM_CPLX(object):
    name = 'SEAM_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX'
    subtables = []


@register_table
class SHEAR(object):
    name = 'SHEAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('F41', '<f8', ()), ('F21', '<f8', ()), ('F12', '<f8', ()), ('F32', '<f8', ()),
             ('F23', '<f8', ()), ('F43', '<f8', ()), ('F34', '<f8', ()), ('F14', '<f8', ()), ('KF1', '<f8', ()),
             ('S12', '<f8', ()), ('KF2', '<f8', ()), ('S23', '<f8', ()), ('KF3', '<f8', ()), ('S34', '<f8', ()),
             ('KF4', '<f8', ()), ('S41', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SHEAR_CPLX(object):
    name = 'SHEAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('F41R', '<f8', ()), ('F21R', '<f8', ()), ('F12R', '<f8', ()), ('F32R', '<f8', ()),
             ('F23R', '<f8', ()), ('F43R', '<f8', ()), ('F34R', '<f8', ()), ('F14R', '<f8', ()), ('KF1R', '<f8', ()),
             ('S12R', '<f8', ()), ('KF2R', '<f8', ()), ('S23R', '<f8', ()), ('KF3R', '<f8', ()), ('S34R', '<f8', ()),
             ('KF4R', '<f8', ()), ('S41R', '<f8', ()), ('F41I', '<f8', ()), ('F21I', '<f8', ()), ('F12I', '<f8', ()),
             ('F32I', '<f8', ()), ('F23I', '<f8', ()), ('F43I', '<f8', ()), ('F34I', '<f8', ()), ('F14I', '<f8', ()),
             ('KF1I', '<f8', ()), ('S12I', '<f8', ()), ('KF2I', '<f8', ()), ('S23I', '<f8', ()), ('KF3I', '<f8', ()),
             ('S34I', '<f8', ()), ('KF4I', '<f8', ()), ('S41I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3(object):
    name = 'TRIA3'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4'
    subtables = []


@register_table
class TRIA3_COMP(object):
    name = 'TRIA3_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_COMP'
    subtables = []


@register_table
class TRIA3_CPLX(object):
    name = 'TRIA3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CPLX'
    subtables = []


@register_table
class TRIA6(object):
    name = 'TRIA6'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('FORCE',
                                                      [('GRID', '<i8', ()), ('MX', '<f8', ()), ('MY', '<f8', ()),
                                                       ('MXY', '<f8', ()), ('BMX', '<f8', ()), ('BMY', '<f8', ()),
                                                       ('BMXY', '<f8', ()), ('TX', '<f8', ()), ('TY', '<f8', ())],
                                                      (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA6_COMP(object):
    name = 'TRIA6_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_COMP'
    subtables = []


@register_table
class TRIA6_CPLX(object):
    name = 'TRIA6_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('FORCE',
                                                      [('GRID', '<i8', ()), ('MXR', '<f8', ()), ('MYR', '<f8', ()),
                                                       ('MXYR', '<f8', ()), ('BMXR', '<f8', ()), ('BMYR', '<f8', ()),
                                                       ('BMXYR', '<f8', ()), ('TXR', '<f8', ()), ('TYR', '<f8', ()),
                                                       ('MXI', '<f8', ()), ('MYI', '<f8', ()), ('MXYI', '<f8', ()),
                                                       ('BMXI', '<f8', ()), ('BMYI', '<f8', ()), ('BMXYI', '<f8', ()),
                                                       ('TXI', '<f8', ()), ('TYI', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIAR(object):
    name = 'TRIAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6'
    subtables = []


@register_table
class TRIAR_CPLX(object):
    name = 'TRIAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6_CPLX'
    subtables = []


@register_table
class TUBE(object):
    name = 'TUBE'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD'
    subtables = []


@register_table
class TUBE_CPLX(object):
    name = 'TUBE_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD_CPLX'
    subtables = []


@register_table
class VISC(object):
    name = 'VISC'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD'
    subtables = []


@register_table
class VISC_CPLX(object):
    name = 'VISC_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/ROD_CPLX'
    subtables = []


@register_table
class WELD(object):
    name = 'WELD'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR'
    subtables = []


@register_table
class WELDC(object):
    name = 'WELDC'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR'
    subtables = []


@register_table
class WELDC_CPLX(object):
    name = 'WELDC_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX'
    subtables = []


@register_table
class WELDP(object):
    name = 'WELDP'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR'
    subtables = []


@register_table
class WELDP_CPLX(object):
    name = 'WELDP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX'
    subtables = []


@register_table
class WELD_CPLX(object):
    name = 'WELD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = [('IDENT', '<i8', ()), ('ELNAME', 'S8', ()), ('ETOTAL', '<f8', ()), ('CVALRES', '<i8', ()),
             ('ESUBT', '<f8', ()), ('ETOTPOS', '<f8', ()), ('ETOTNEG', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class KINETIC_DMIG(object):
    name = 'KINETIC_DMIG'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ENERGY/STRAIN_DMIG'
    subtables = []


@register_table
class KINETIC_ELEM(object):
    name = 'KINETIC_ELEM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ENERGY/STRAIN_ELEM'
    subtables = []


@register_table
class LOSS_DMIG(object):
    name = 'LOSS_DMIG'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ENERGY/STRAIN_DMIG'
    subtables = []


@register_table
class LOSS_ELEM(object):
    name = 'LOSS_ELEM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/ENERGY/STRAIN_ELEM'
    subtables = []


@register_table
class STRAIN_DMIG(object):
    name = 'STRAIN_DMIG'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = [('NAME', 'S8', ()), ('ENERGY', '<f8', ()), ('PCT', '<f8', ()), ('DEN', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class STRAIN_ELEM(object):
    name = 'STRAIN_ELEM'
    path = '/NASTRAN/RESULT/ELEMENTAL/ENERGY'
    dtype = [('ID', '<i8', ()), ('ENERGY', '<f8', ()), ('PCT', '<f8', ()), ('DEN', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DYNAMIC(object):
    name = 'DYNAMIC'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY/MASS'
    subtables = []


@register_table
class MASS(object):
    name = 'MASS'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = [('EID', '<i8', ()), ('ELTYPE', 'S8', ()), ('RESPDOF', '<i8', ()), ('SENSR', '<f8', ()),
             ('SENSI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SQ_DYNAMIC(object):
    name = 'SQ_DYNAMIC'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY/MASS'
    subtables = []


@register_table
class SQ_MASS(object):
    name = 'SQ_MASS'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY/MASS'
    subtables = []


@register_table
class SQ_STIFFNESS(object):
    name = 'SQ_STIFFNESS'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY/MASS'
    subtables = []


@register_table
class STIFFNESS(object):
    name = 'STIFFNESS'
    path = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/SENSITIVITY/MASS'
    subtables = []


@register_table
class AXIF2(object):
    name = 'AXIF2'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF2'
    subtables = []


@register_table
class AXIF2_CPLX(object):
    name = 'AXIF2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF2_CPLX'
    subtables = []


@register_table
class AXIF3(object):
    name = 'AXIF3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF3'
    subtables = []


@register_table
class AXIF3_CPLX(object):
    name = 'AXIF3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF3_CPLX'
    subtables = []


@register_table
class AXIF4(object):
    name = 'AXIF4'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF4'
    subtables = []


@register_table
class AXIF4_CPLX(object):
    name = 'AXIF4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/AXIF4_CPLX'
    subtables = []


@register_table
class BAR(object):
    name = 'BAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BAR'
    subtables = []


@register_table
class BARS(object):
    name = 'BARS'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BARS'
    subtables = []


@register_table
class BARS_CPLX(object):
    name = 'BARS_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BARS_CPLX'
    subtables = []


@register_table
class BARS_RR(object):
    name = 'BARS_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BARS_RR'
    subtables = []


@register_table
class BAR_CPLX(object):
    name = 'BAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BAR_CPLX'
    subtables = []


@register_table
class BAR_RR(object):
    name = 'BAR_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BAR_RR'
    subtables = []


@register_table
class BEAM(object):
    name = 'BEAM'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM'
    subtables = []


@register_table
class BEAM3(object):
    name = 'BEAM3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM3'
    subtables = []


@register_table
class BEAM3_CPLX(object):
    name = 'BEAM3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM3_CPLX'
    subtables = []


@register_table
class BEAM3_RR(object):
    name = 'BEAM3_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM3_RR'
    subtables = []


@register_table
class BEAM_CPLX(object):
    name = 'BEAM_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM_CPLX'
    subtables = []


@register_table
class BEAM_RR(object):
    name = 'BEAM_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEAM_RR'
    subtables = []


@register_table
class BEND(object):
    name = 'BEND'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEND'
    subtables = []


@register_table
class BEND_CPLX(object):
    name = 'BEND_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEND_CPLX'
    subtables = []


@register_table
class BEND_RR(object):
    name = 'BEND_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BEND_RR'
    subtables = []


@register_table
class BUSH(object):
    name = 'BUSH'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BUSH'
    subtables = []


@register_table
class BUSH1D(object):
    name = 'BUSH1D'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BUSH1D'
    subtables = []


@register_table
class BUSH1D_CPLX(object):
    name = 'BUSH1D_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BUSH1D_CPLX'
    subtables = []


@register_table
class BUSH1D_RR(object):
    name = 'BUSH1D_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BUSH1D_RR'
    subtables = []


@register_table
class BUSH_CPLX(object):
    name = 'BUSH_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/BUSH_CPLX'
    subtables = []


@register_table
class CONE(object):
    name = 'CONE'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/CONE'
    subtables = []


@register_table
class CONROD(object):
    name = 'CONROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/CONROD'
    subtables = []


@register_table
class CONROD_CPLX(object):
    name = 'CONROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/CONROD_CPLX'
    subtables = []


@register_table
class CONROD_RR(object):
    name = 'CONROD_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/CONROD_RR'
    subtables = []


@register_table
class ELAS1(object):
    name = 'ELAS1'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1'
    subtables = []


@register_table
class ELAS1_CPLX(object):
    name = 'ELAS1_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS2(object):
    name = 'ELAS2'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1'
    subtables = []


@register_table
class ELAS2_CPLX(object):
    name = 'ELAS2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS3(object):
    name = 'ELAS3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1'
    subtables = []


@register_table
class ELAS3_CPLX(object):
    name = 'ELAS3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_CPLX'
    subtables = []


@register_table
class BEAM3(object):
    name = 'BEAM3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN/ELEM_COMP'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/BEAM3'
    subtables = []


@register_table
class FAST(object):
    name = 'FAST'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/FAST'
    subtables = []


@register_table
class FAST_CPLX(object):
    name = 'FAST_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/FAST_CPLX'
    subtables = []


@register_table
class GAP(object):
    name = 'GAP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/GAP'
    subtables = []


@register_table
class HEXA(object):
    name = 'HEXA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/HEXA'
    subtables = []


@register_table
class HEXA_CPLX(object):
    name = 'HEXA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/HEXA_CPLX'
    subtables = []


@register_table
class PENTA(object):
    name = 'PENTA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/PENTA'
    subtables = []


@register_table
class PENTA_CPLX(object):
    name = 'PENTA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/PENTA_CPLX'
    subtables = []


@register_table
class QUAD4(object):
    name = 'QUAD4'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4'
    subtables = []


@register_table
class QUAD4_COMP(object):
    name = 'QUAD4_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class QUAD4_COMP_CPLX(object):
    name = 'QUAD4_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class QUAD4_CPLX(object):
    name = 'QUAD4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_CPLX'
    subtables = []


@register_table
class QUAD4_FD(object):
    name = 'QUAD4_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FD'
    subtables = []


@register_table
class QUAD4_FD_CPLX(object):
    name = 'QUAD4_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FD_CPLX'
    subtables = []


@register_table
class QUAD8(object):
    name = 'QUAD8'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8'
    subtables = []


@register_table
class QUAD8_COMP(object):
    name = 'QUAD8_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_COMP'
    subtables = []


@register_table
class QUAD8_COMP_CPLX(object):
    name = 'QUAD8_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_COMP_CPLX'
    subtables = []


@register_table
class QUAD8_CPLX(object):
    name = 'QUAD8_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_CPLX'
    subtables = []


@register_table
class QUAD8_FD(object):
    name = 'QUAD8_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_FD'
    subtables = []


@register_table
class QUAD8_FD_CPLX(object):
    name = 'QUAD8_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_FD_CPLX'
    subtables = []


@register_table
class QUADR(object):
    name = 'QUADR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR'
    subtables = []


@register_table
class QUADR_COMP(object):
    name = 'QUADR_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR_COMP'
    subtables = []


@register_table
class QUADR_COMP_CPLX(object):
    name = 'QUADR_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR_COMP_CPLX'
    subtables = []


@register_table
class QUADR_CPLX(object):
    name = 'QUADR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR_CPLX'
    subtables = []


@register_table
class QUADR_FD(object):
    name = 'QUADR_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR_FD'
    subtables = []


@register_table
class QUADR_FD_CPLX(object):
    name = 'QUADR_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADR_FD_CPLX'
    subtables = []


@register_table
class QUADX4_FD(object):
    name = 'QUADX4_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADX4_FD'
    subtables = []


@register_table
class QUADX4_FD_CPLX(object):
    name = 'QUADX4_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADX4_FD_CPLX'
    subtables = []


@register_table
class QUADX8_FD(object):
    name = 'QUADX8_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADX8_FD'
    subtables = []


@register_table
class QUADX8_FD_CPLX(object):
    name = 'QUADX8_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUADX8_FD_CPLX'
    subtables = []


@register_table
class QUAD_CN(object):
    name = 'QUAD_CN'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN'
    subtables = []


@register_table
class QUAD_CN_CPLX(object):
    name = 'QUAD_CN_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN_CPLX'
    subtables = []


@register_table
class ROD(object):
    name = 'ROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ROD'
    subtables = []


@register_table
class ROD_CPLX(object):
    name = 'ROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ROD_CPLX'
    subtables = []


@register_table
class ROD_RR(object):
    name = 'ROD_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ROD_RR'
    subtables = []


@register_table
class SEAMP(object):
    name = 'SEAMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/SEAMP'
    subtables = []


@register_table
class SEAMP_CPLX(object):
    name = 'SEAMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/SEAMP_CPLX'
    subtables = []


@register_table
class SHEAR(object):
    name = 'SHEAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/SHEAR'
    subtables = []


@register_table
class SHEAR_CPLX(object):
    name = 'SHEAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/SHEAR_CPLX'
    subtables = []


@register_table
class SHEAR_RR(object):
    name = 'SHEAR_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/SHEAR_RR'
    subtables = []


@register_table
class TETRA(object):
    name = 'TETRA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TETRA'
    subtables = []


@register_table
class TETRA_CPLX(object):
    name = 'TETRA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TETRA_CPLX'
    subtables = []


@register_table
class TRIA3(object):
    name = 'TRIA3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3'
    subtables = []


@register_table
class TRIA3_COMP(object):
    name = 'TRIA3_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_COMP'
    subtables = []


@register_table
class TRIA3_COMP_CPLX(object):
    name = 'TRIA3_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_COMP_CPLX'
    subtables = []


@register_table
class TRIA3_CPLX(object):
    name = 'TRIA3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_CPLX'
    subtables = []


@register_table
class TRIA3_FD(object):
    name = 'TRIA3_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_FD'
    subtables = []


@register_table
class TRIA3_FD_CPLX(object):
    name = 'TRIA3_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_FD_CPLX'
    subtables = []


@register_table
class TRIA6(object):
    name = 'TRIA6'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6'
    subtables = []


@register_table
class TRIA6_COMP(object):
    name = 'TRIA6_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_COMP'
    subtables = []


@register_table
class TRIA6_COMP_CPLX(object):
    name = 'TRIA6_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_COMP_CPLX'
    subtables = []


@register_table
class TRIA6_CPLX(object):
    name = 'TRIA6_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_CPLX'
    subtables = []


@register_table
class TRIA6_FD(object):
    name = 'TRIA6_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_FD'
    subtables = []


@register_table
class TRIA6_FD_CPLX(object):
    name = 'TRIA6_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_FD_CPLX'
    subtables = []


@register_table
class TRIAR(object):
    name = 'TRIAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR'
    subtables = []


@register_table
class TRIAR_1FD(object):
    name = 'TRIAR_1FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_1FD'
    subtables = []


@register_table
class TRIAR_1FD_CPLX(object):
    name = 'TRIAR_1FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_1FD_CPLX'
    subtables = []


@register_table
class TRIAR_4FD(object):
    name = 'TRIAR_4FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_4FD'
    subtables = []


@register_table
class TRIAR_4FD_CPLX(object):
    name = 'TRIAR_4FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_4FD_CPLX'
    subtables = []


@register_table
class TRIAR_COMP(object):
    name = 'TRIAR_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_COMP'
    subtables = []


@register_table
class TRIAR_COMP_CPLX(object):
    name = 'TRIAR_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_COMP_CPLX'
    subtables = []


@register_table
class TRIAR_CPLX(object):
    name = 'TRIAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAR_CPLX'
    subtables = []


@register_table
class TRIAX3_FD(object):
    name = 'TRIAX3_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX3_FD'
    subtables = []


@register_table
class TRIAX3_FD_CPLX(object):
    name = 'TRIAX3_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX3_FD_CPLX'
    subtables = []


@register_table
class TRIAX6(object):
    name = 'TRIAX6'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX6'
    subtables = []


@register_table
class TRIAX6_CPLX(object):
    name = 'TRIAX6_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX6_CPLX'
    subtables = []


@register_table
class TRIAX6_FD(object):
    name = 'TRIAX6_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX6_FD'
    subtables = []


@register_table
class TRIAX6_FD_CPLX(object):
    name = 'TRIAX6_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIAX6_FD_CPLX'
    subtables = []


@register_table
class TUBE(object):
    name = 'TUBE'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TUBE'
    subtables = []


@register_table
class TUBE_CPLX(object):
    name = 'TUBE_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TUBE_CPLX'
    subtables = []


@register_table
class TUBE_RR(object):
    name = 'TUBE_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TUBE_RR'
    subtables = []


@register_table
class WELD(object):
    name = 'WELD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD'
    subtables = []


@register_table
class WELD_CPLX(object):
    name = 'WELD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD_CPLX'
    subtables = []


@register_table
class WELDC(object):
    name = 'WELDC'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD'
    subtables = []


@register_table
class WELDC_CPLX(object):
    name = 'WELDC_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD_CPLX'
    subtables = []


@register_table
class WELDP(object):
    name = 'WELDP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD'
    subtables = []


@register_table
class WELDP_CPLX(object):
    name = 'WELDP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRAIN'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD_CPLX'
    subtables = []


@register_table
class AXIF2(object):
    name = 'AXIF2'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RA', '<f8', ()), ('AA', '<f8', ()), ('TE', '<f8', ()), ('CE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIF2_CPLX(object):
    name = 'AXIF2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RAR', '<f8', ()), ('AAR', '<f8', ()), ('TER', '<f8', ()), ('CER', '<f8', ()),
             ('RAI', '<f8', ()), ('AAI', '<f8', ()), ('TEI', '<f8', ()), ('CEI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIF3(object):
    name = 'AXIF3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RC', '<f8', ()), ('CC', '<f8', ()), ('AC', '<f8', ()), ('TE1', '<f8', ()),
             ('CE1', '<f8', ()), ('TE2', '<f8', ()), ('CE2', '<f8', ()), ('TE3', '<f8', ()), ('CE3', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIF3_CPLX(object):
    name = 'AXIF3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RCR', '<f8', ()), ('CCR', '<f8', ()), ('ACR', '<f8', ()), ('TE1R', '<f8', ()),
             ('CE1R', '<f8', ()), ('TE2R', '<f8', ()), ('CE2R', '<f8', ()), ('TE3R', '<f8', ()), ('CE3R', '<f8', ()),
             ('RCI', '<f8', ()), ('CCI', '<f8', ()), ('ACI', '<f8', ()), ('TE1I', '<f8', ()), ('CE1I', '<f8', ()),
             ('TE2I', '<f8', ()), ('CE2I', '<f8', ()), ('TE3I', '<f8', ()), ('CE3I', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIF4(object):
    name = 'AXIF4'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RC', '<f8', ()), ('CC', '<f8', ()), ('AC', '<f8', ()), ('TE1', '<f8', ()),
             ('CE1', '<f8', ()), ('TE2', '<f8', ()), ('CE2', '<f8', ()), ('TE3', '<f8', ()), ('CE3', '<f8', ()),
             ('TE4', '<f8', ()), ('CE4', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXIF4_CPLX(object):
    name = 'AXIF4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RCR', '<f8', ()), ('CCR', '<f8', ()), ('ACR', '<f8', ()), ('TE1R', '<f8', ()),
             ('CE1R', '<f8', ()), ('TE2R', '<f8', ()), ('CE2R', '<f8', ()), ('TE3R', '<f8', ()), ('CE3R', '<f8', ()),
             ('TE4R', '<f8', ()), ('CE4R', '<f8', ()), ('RCI', '<f8', ()), ('CCI', '<f8', ()), ('ACI', '<f8', ()),
             ('TE1I', '<f8', ()), ('CE1I', '<f8', ()), ('TE2I', '<f8', ()), ('CE2I', '<f8', ()), ('TE3I', '<f8', ()),
             ('CE3I', '<f8', ()), ('TE4I', '<f8', ()), ('CE4I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AXISYM(object):
    name = 'AXISYM'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FD1', '<f8', ()), ('SXX1', '<f8', ()), ('SYY1', '<f8', ()), ('SXY1', '<f8', ()),
             ('EXX1', '<f8', ()), ('EYY1', '<f8', ()), ('EXY1', '<f8', ()), ('FD2', '<f8', ()), ('SXX2', '<f8', ()),
             ('SYY2', '<f8', ()), ('SXY2', '<f8', ()), ('EXX2', '<f8', ()), ('EYY2', '<f8', ()), ('EXY2', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR(object):
    name = 'BAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('X1A', '<f8', ()), ('X2A', '<f8', ()), ('X3A', '<f8', ()), ('X4A', '<f8', ()),
             ('AX', '<f8', ()), ('MAXA', '<f8', ()), ('MINA', '<f8', ()), ('MST', '<f8', ()), ('X1B', '<f8', ()),
             ('X2B', '<f8', ()), ('X3B', '<f8', ()), ('X4B', '<f8', ()), ('MAXB', '<f8', ()), ('MINB', '<f8', ()),
             ('MSC', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BARS(object):
    name = 'BARS'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SD', '<f8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()), ('XE', '<f8', ()),
             ('XF', '<f8', ()), ('AX', '<f8', ()), ('MAX', '<f8', ()), ('MIN', '<f8', ()), ('MS', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BARS_CPLX(object):
    name = 'BARS_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SD', '<f8', ()), ('XCR', '<f8', ()), ('XDR', '<f8', ()), ('XER', '<f8', ()),
             ('XFR', '<f8', ()), ('AXR', '<f8', ()), ('MAXR', '<f8', ()), ('MINR', '<f8', ()), ('XCI', '<f8', ()),
             ('XDI', '<f8', ()), ('XEI', '<f8', ()), ('XFI', '<f8', ()), ('AXI', '<f8', ()), ('MAXI', '<f8', ()),
             ('MINI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BARS_RR(object):
    name = 'BARS_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SD', '<f8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()), ('XE', '<f8', ()),
             ('XF', '<f8', ()), ('AX', '<f8', ()), ('MAX', '<f8', ()), ('MIN', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR_CPLX(object):
    name = 'BAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('X1AR', '<f8', ()), ('X2AR', '<f8', ()), ('X3AR', '<f8', ()), ('X4AR', '<f8', ()),
             ('AXR', '<f8', ()), ('X1AI', '<f8', ()), ('X2AI', '<f8', ()), ('X3AI', '<f8', ()), ('X4AI', '<f8', ()),
             ('AXI', '<f8', ()), ('X1BR', '<f8', ()), ('X2BR', '<f8', ()), ('X3BR', '<f8', ()), ('X4BR', '<f8', ()),
             ('X1BI', '<f8', ()), ('X2BI', '<f8', ()), ('X3BI', '<f8', ()), ('X4BI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR_NL(object):
    name = 'BAR_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('GRIDA', '<i8', ()), ('NSXCA', '<f8', ()), ('NSECA', '<f8', ()), ('TECA', '<f8', ()),
             ('EPECA', '<f8', ()), ('ECECA', '<f8', ()), ('NSXDA', '<f8', ()), ('NSEDA', '<f8', ()),
             ('TEDA', '<f8', ()), ('EPEDA', '<f8', ()), ('ECEDA', '<f8', ()), ('NSXEA', '<f8', ()),
             ('NSEEA', '<f8', ()), ('TEEA', '<f8', ()), ('EPEEA', '<f8', ()), ('ECEEA', '<f8', ()),
             ('NSXFA', '<f8', ()), ('NSEFA', '<f8', ()), ('TEFA', '<f8', ()), ('EPEFA', '<f8', ()),
             ('ECEFA', '<f8', ()), ('GRIDB', '<i8', ()), ('NSXCB', '<f8', ()), ('NSECB', '<f8', ()),
             ('TECB', '<f8', ()), ('EPECB', '<f8', ()), ('ECECB', '<f8', ()), ('NSXDB', '<f8', ()),
             ('NSEDB', '<f8', ()), ('TEDB', '<f8', ()), ('EPEDB', '<f8', ()), ('ECEDB', '<f8', ()),
             ('NSXEB', '<f8', ()), ('NSEEB', '<f8', ()), ('TEEB', '<f8', ()), ('EPEEB', '<f8', ()),
             ('ECEEB', '<f8', ()), ('NSXFB', '<f8', ()), ('NSEFB', '<f8', ()), ('TEFB', '<f8', ()),
             ('EPEFB', '<f8', ()), ('ECEFB', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BAR_RR(object):
    name = 'BAR_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('X1A', '<f8', ()), ('X2A', '<f8', ()), ('X3A', '<f8', ()), ('X4A', '<f8', ()),
             ('AX', '<f8', ()), ('X1B', '<f8', ()), ('X2B', '<f8', ()), ('X3B', '<f8', ()), ('X4B', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM(object):
    name = 'BEAM'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('SD', '<f8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()),
                                         ('XE', '<f8', ()), ('XF', '<f8', ()), ('MAX', '<f8', ()), ('MIN', '<f8', ()),
                                         ('MST', '<f8', ()), ('MSC', '<f8', ())], (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM3(object):
    name = 'BEAM3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()), ('XE', '<f8', ()),
                                         ('XF', '<f8', ()), ('MAX', '<f8', ()), ('MIN', '<f8', ()), ('MST', '<f8', ()),
                                         ('MSC', '<f8', ()), ('YC', '<f8', ()), ('YD', '<f8', ()), ('YE', '<f8', ()),
                                         ('YF', '<f8', ()), ('YMAX', '<f8', ()), ('YMIN', '<f8', ()), ('ZC', '<f8', ()),
                                         ('ZD', '<f8', ()), ('ZE', '<f8', ()), ('ZF', '<f8', ()), ('ZMAX', '<f8', ()),
                                         ('ZMIN', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM3_CPLX(object):
    name = 'BEAM3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('XCR', '<f8', ()), ('XDR', '<f8', ()), ('XER', '<f8', ()),
                                   ('XFR', '<f8', ()), ('XCI', '<f8', ()), ('XDI', '<f8', ()), ('XEI', '<f8', ()),
                                   ('XFI', '<f8', ()), ('YCR', '<f8', ()), ('YDR', '<f8', ()), ('YER', '<f8', ()),
                                   ('YFR', '<f8', ()), ('YCI', '<f8', ()), ('YDI', '<f8', ()), ('YEI', '<f8', ()),
                                   ('YFI', '<f8', ()), ('ZCR', '<f8', ()), ('ZDR', '<f8', ()), ('ZER', '<f8', ()),
                                   ('ZFR', '<f8', ()), ('ZCI', '<f8', ()), ('ZDI', '<f8', ()), ('ZEI', '<f8', ()),
                                   ('ZFI', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM3_RR(object):
    name = 'BEAM3_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()), ('XE', '<f8', ()),
                                         ('XF', '<f8', ()), ('YC', '<f8', ()), ('YD', '<f8', ()), ('YE', '<f8', ()),
                                         ('YF', '<f8', ()), ('ZC', '<f8', ()), ('ZD', '<f8', ()), ('ZE', '<f8', ()),
                                         ('ZF', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_CPLX(object):
    name = 'BEAM_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('SD', '<f8', ()), ('XCR', '<f8', ()), ('XDR', '<f8', ()),
                                         ('XER', '<f8', ()), ('XFR', '<f8', ()), ('XCI', '<f8', ()), ('XDI', '<f8', ()),
                                         ('XEI', '<f8', ()), ('XFI', '<f8', ())], (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_NL(object):
    name = 'BEAM_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('GRIDA', '<i8', ()), ('NSXCA', '<f8', ()), ('NSECA', '<f8', ()), ('TECA', '<f8', ()),
             ('EPECA', '<f8', ()), ('ECECA', '<f8', ()), ('NSXDA', '<f8', ()), ('NSEDA', '<f8', ()),
             ('TEDA', '<f8', ()), ('EPEDA', '<f8', ()), ('ECEDA', '<f8', ()), ('NSXEA', '<f8', ()),
             ('NSEEA', '<f8', ()), ('TEEA', '<f8', ()), ('EPEEA', '<f8', ()), ('ECEEA', '<f8', ()),
             ('NSXFA', '<f8', ()), ('NSEFA', '<f8', ()), ('TEFA', '<f8', ()), ('EPEFA', '<f8', ()),
             ('ECEFA', '<f8', ()), ('GRIDB', '<i8', ()), ('NSXCB', '<f8', ()), ('NSECB', '<f8', ()),
             ('TECB', '<f8', ()), ('EPECB', '<f8', ()), ('ECECB', '<f8', ()), ('NSXDB', '<f8', ()),
             ('NSEDB', '<f8', ()), ('TEDB', '<f8', ()), ('EPEDB', '<f8', ()), ('ECEDB', '<f8', ()),
             ('NSXEB', '<f8', ()), ('NSEEB', '<f8', ()), ('TEEB', '<f8', ()), ('EPEEB', '<f8', ()),
             ('ECEEB', '<f8', ()), ('NSXFB', '<f8', ()), ('NSEFB', '<f8', ()), ('TEFB', '<f8', ()),
             ('EPEFB', '<f8', ()), ('ECEFB', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_RR(object):
    name = 'BEAM_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('SD', '<f8', ()), ('XC', '<f8', ()), ('XD', '<f8', ()),
                                         ('XE', '<f8', ()), ('XF', '<f8', ())], (11,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEND(object):
    name = 'BEND'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('CA', '<f8', ()), ('C', '<f8', ()), ('D', '<f8', ()),
                                         ('E', '<f8', ()), ('F', '<f8', ()), ('MAX', '<f8', ()), ('MIN', '<f8', ()),
                                         ('MST', '<f8', ()), ('MSC', '<f8', ())], (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEND_CPLX(object):
    name = 'BEND_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('CA', '<f8', ()), ('CR', '<f8', ()), ('DR', '<f8', ()),
                                         ('ER', '<f8', ()), ('FR', '<f8', ()), ('CI', '<f8', ()), ('DI', '<f8', ()),
                                         ('EI', '<f8', ()), ('FI', '<f8', ())], (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEND_RR(object):
    name = 'BEND_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('CA', '<f8', ()), ('C', '<f8', ()), ('D', '<f8', ()),
                                         ('E', '<f8', ()), ('F', '<f8', ())], (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH(object):
    name = 'BUSH'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TX', '<f8', ()), ('TY', '<f8', ()), ('TZ', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH1D(object):
    name = 'BUSH1D'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FE', '<f8', ()), ('UE', '<f8', ()), ('VE', '<f8', ()), ('AS', '<f8', ()),
             ('AE', '<f8', ()), ('EP', '<f8', ()), ('FAIL', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH1D_CPLX(object):
    name = 'BUSH1D_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FER', '<f8', ()), ('UER', '<f8', ()), ('ASR', '<f8', ()), ('AER', '<f8', ()),
             ('FEI', '<f8', ()), ('UEI', '<f8', ()), ('ASI', '<f8', ()), ('AEI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH1D_RR(object):
    name = 'BUSH1D_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FE', '<f8', ()), ('UE', '<f8', ()), ('AS', '<f8', ()), ('AE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH_CPLX(object):
    name = 'BUSH_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TXR', '<f8', ()), ('TYR', '<f8', ()), ('TZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('TXI', '<f8', ()), ('TYI', '<f8', ()), ('TZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BUSH_NL(object):
    name = 'BUSH_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FX', '<f8', ()), ('FY', '<f8', ()), ('FZ', '<f8', ()), ('STX', '<f8', ()),
             ('STY', '<f8', ()), ('STZ', '<f8', ()), ('ETX', '<f8', ()), ('ETY', '<f8', ()), ('ETZ', '<f8', ()),
             ('MX', '<f8', ()), ('MY', '<f8', ()), ('MZ', '<f8', ()), ('SRX', '<f8', ()), ('SRY', '<f8', ()),
             ('SRZ', '<f8', ()), ('ERX', '<f8', ()), ('ERY', '<f8', ()), ('ERZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONE(object):
    name = 'CONE'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('HOPA', '<f8', ()), ('FD1', '<f8', ()), ('U1', '<f8', ()), ('V1', '<f8', ()),
             ('T1', '<f8', ()), ('A1', '<f8', ()), ('MJRP1', '<f8', ()), ('MNRP1', '<f8', ()), ('MAX1', '<f8', ()),
             ('FD2', '<f8', ()), ('U2', '<f8', ()), ('V2', '<f8', ()), ('T2', '<f8', ()), ('A2', '<f8', ()),
             ('MJRP2', '<f8', ()), ('MNRP2', '<f8', ()), ('MAX2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD(object):
    name = 'CONROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('A', '<f8', ()), ('MSA', '<f8', ()), ('T', '<f8', ()), ('MST', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD_CPLX(object):
    name = 'CONROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AR', '<f8', ()), ('AI', '<f8', ()), ('TR', '<f8', ()), ('TI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD_NL(object):
    name = 'CONROD_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('SE', '<f8', ()), ('TE', '<f8', ()), ('EPS', '<f8', ()),
             ('ECS', '<f8', ()), ('LTS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONROD_RR(object):
    name = 'CONROD_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('A', '<f8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS1(object):
    name = 'ELAS1'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('S', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS1_CPLX(object):
    name = 'ELAS1_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SR', '<f8', ()), ('SI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS1_NL(object):
    name = 'ELAS1_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('F', '<f8', ()), ('S', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELAS2(object):
    name = 'ELAS2'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1'
    subtables = []


@register_table
class ELAS2_CPLX(object):
    name = 'ELAS2_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS3(object):
    name = 'ELAS3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1'
    subtables = []


@register_table
class ELAS3_CPLX(object):
    name = 'ELAS3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_CPLX'
    subtables = []


@register_table
class ELAS3_NL(object):
    name = 'ELAS3_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELAS1_NL'
    subtables = []


@register_table
class BEAM3(object):
    name = 'BEAM3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GID', '<i8', ()), ('PLYID', '<i8', ()), ('NFT', '<i8', ()), ('FT1', 'S8', ()),
                                   ('FT2', 'S8', ()), ('FT3', 'S8', ()), ('SXX', '<f8', ()), ('SYY', '<f8', ()),
                                   ('SZZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()), ('SXZ', '<f8', ()),
                                   ('FP1', '<f8', ()), ('FP2', '<f8', ()), ('FP3', '<f8', ()), ('SR1', '<f8', ()),
                                   ('SR2', '<f8', ()), ('SR3', '<f8', ()), ('FAIL', '<i8', ())], (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM'
    dtype = [('SS',
              [('INTPID', '<i8', ()), ('SXX', '<f8', ()), ('SYY', '<f8', ()), ('SZZ', '<f8', ()), ('SXY', '<f8', ()),
               ('SYZ', '<f8', ()), ('SXZ', '<f8', ()), ('EXX', '<f8', ()), ('EYY', '<f8', ()), ('EZZ', '<f8', ()),
               ('EXY', '<f8', ()), ('EYZ', '<f8', ()), ('EXZ', '<f8', ()), ('SLAMN1', '<f8', ()), ('SLAMN2', '<f8', ()),
               ('SLAMN3', '<f8', ()), ('TLAMN1', '<f8', ()), ('TLAMN2', '<f8', ()), ('TLAMN3', '<f8', ()),
               ('FP1', '<f8', ()), ('FP2', '<f8', ()), ('FP3', '<f8', ()), ('SR1', '<f8', ()), ('SR2', '<f8', ()),
               ('SR3', '<f8', ()), ('FM1', '<i8', ()), ('FM2', '<i8', ()), ('FM3', '<i8', ()), ('FBINDX', '<f8', ()),
               ('DAM1', '<f8', ()), ('DAM2', '<f8', ()), ('DAM3', '<f8', ()), ('DAMM', '<i8', ()),
               ('CREDENS', '<f8', ()), ('DAMM2', '<i8', ())], ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('NFT', '<i8', ()), ('FT11', 'S4', ()), ('FT12', 'S4', ()),
             ('FT21', 'S4', ()), ('FT22', 'S4', ()), ('FT31', 'S4', ()), ('FT32', 'S4', ()), ('DATA_POS', '<i8', ()),
             ('DATA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM/DATA']


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM_CPLX'
    dtype = [('SS', [('INTPID', '<i8', ()), ('SXXR', '<f8', ()), ('SYYR', '<f8', ()), ('SZZR', '<f8', ()),
                     ('SXYR', '<f8', ()), ('SYZR', '<f8', ()), ('SXZR', '<f8', ()), ('EXXR', '<f8', ()),
                     ('EYYR', '<f8', ()), ('EZZR', '<f8', ()), ('EXYR', '<f8', ()), ('EYZR', '<f8', ()),
                     ('EXZR', '<f8', ()), ('SXXI', '<f8', ()), ('SYYI', '<f8', ()), ('SZZI', '<f8', ()),
                     ('SXYI', '<f8', ()), ('SYZI', '<f8', ()), ('SXZI', '<f8', ()), ('EXXI', '<f8', ()),
                     ('EYYI', '<f8', ()), ('EZZI', '<f8', ()), ('EXYI', '<f8', ()), ('EYZI', '<f8', ()),
                     ('EXZI', '<f8', ()), ('SLAMN1R', '<f8', ()), ('SLAMN2R', '<f8', ()), ('SLAMN3R', '<f8', ()),
                     ('TLAMN1R', '<f8', ()), ('TLAMN2R', '<f8', ()), ('TLAMN3R', '<f8', ()), ('SLAMN1I', '<f8', ()),
                     ('SLAMN2I', '<f8', ()), ('SLAMN3I', '<f8', ()), ('TLAMN1I', '<f8', ()), ('TLAMN2I', '<f8', ()),
                     ('TLAMN3I', '<f8', ()), ('FP1', '<f8', ()), ('FP2', '<f8', ()), ('FP3', '<f8', ()),
                     ('SR1', '<f8', ()), ('SR2', '<f8', ()), ('SR3', '<f8', ()), ('FM1', '<i8', ()), ('FM2', '<i8', ()),
                     ('FM3', '<i8', ()), ('FBINDX', '<f8', ())], ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM_CPLX'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('NFT', '<i8', ()), ('FT11', 'S4', ()), ('FT12', 'S4', ()),
             ('FT21', 'S4', ()), ('FT22', 'S4', ()), ('FT31', 'S4', ()), ('FT32', 'S4', ()), ('DATA_POS', '<i8', ()),
             ('DATA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/ELEM_CPLX/DATA']


@register_table
class MAX_MIN(object):
    name = 'MAX_MIN'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP'
    dtype = [('EID', '<i8', ()), ('ELTYPE', '<i8', ()), ('FMAX', '<f8', ()), ('PLYFMAX', '<i8', ()),
             ('SRMIN', '<f8', ()), ('PLYSRMIN', '<i8', ()), ('FAILFLG', 'S4', ()), ('DAMAGE', '<i8', ()),
             ('DMAX', '<f8', ()), ('PLYDMAX', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/GASKET'
    dtype = [('SS', [('INTPID', '<i8', ()), ('GPRES', '<f8', ()), ('GCLOS', '<f8', ()), ('GPCLOS', '<f8', ())], ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/GASKET'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('DATA_POS', '<i8', ()), ('DATA_LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/GASKET/DATA']


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/THERM'
    dtype = [('SS', [('INTPID', '<i8', ()), ('GRADX', '<f8', ()), ('GRADY', '<f8', ()), ('GRADZ', '<f8', ()),
                     ('FLUXX', '<f8', ()), ('FLUXY', '<f8', ()), ('FLUXZ', '<f8', ()), ('TEMP', '<f8', ())], ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/THERM'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('THICK', '<f8', ()), ('DATA_POS', '<i8', ()),
             ('DATA_LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/ELEMENTAL/STRESS/ELEM_COMP/THERM/DATA']


@register_table
class EXTREME_FIBRE(object):
    name = 'EXTREME_FIBRE'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('ID', '<i8', ()), ('FT', '<i8', ()), ('PLY1', '<i8', ()), ('SX1', '<f8', ()), ('SY1', '<f8', ()),
             ('T1', '<f8', ()), ('SL11', '<f8', ()), ('SL21', '<f8', ()), ('FP1', '<f8', ()), ('FB1', '<f8', ()),
             ('EX1', '<f8', ()), ('EY1', '<f8', ()), ('ET1', '<f8', ()), ('EL11', '<f8', ()), ('EL21', '<f8', ()),
             ('SRP1', '<f8', ()), ('SRB1', '<f8', ()), ('PLY2', '<i8', ()), ('SX2', '<f8', ()), ('SY2', '<f8', ()),
             ('T2', '<f8', ()), ('SL12', '<f8', ()), ('SL22', '<f8', ()), ('FP2', '<f8', ()), ('FB2', '<f8', ()),
             ('EX2', '<f8', ()), ('EY2', '<f8', ()), ('ET2', '<f8', ()), ('EL12', '<f8', ()), ('EL22', '<f8', ()),
             ('SRP2', '<f8', ()), ('SRB2', '<f8', ()), ('ELTYPE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EXTREME_FIBRE_CPLX(object):
    name = 'EXTREME_FIBRE_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('ID', '<i8', ()), ('FT', '<i8', ()), ('PLY1', '<i8', ()), ('SX1R', '<f8', ()), ('SX1I', '<f8', ()),
             ('SY1R', '<f8', ()), ('SY1I', '<f8', ()), ('T1R', '<f8', ()), ('T1I', '<f8', ()), ('SL11R', '<f8', ()),
             ('SL11I', '<f8', ()), ('SL21R', '<f8', ()), ('SL21I', '<f8', ()), ('FP1R', '<f8', ()), ('FP1I', '<f8', ()),
             ('FB1R', '<f8', ()), ('FB1I', '<f8', ()), ('EX1R', '<f8', ()), ('EX1I', '<f8', ()), ('EY1R', '<f8', ()),
             ('EY1I', '<f8', ()), ('ET1R', '<f8', ()), ('ET1I', '<f8', ()), ('EL11R', '<f8', ()), ('EL11I', '<f8', ()),
             ('EL21R', '<f8', ()), ('EL21I', '<f8', ()), ('SRP1R', '<f8', ()), ('SRP1I', '<f8', ()),
             ('SRB1R', '<f8', ()), ('SRB1I', '<f8', ()), ('PLY2', '<i8', ()), ('SX2R', '<f8', ()), ('SX2I', '<f8', ()),
             ('SY2R', '<f8', ()), ('SY2I', '<f8', ()), ('T2R', '<f8', ()), ('T2I', '<f8', ()), ('SL12R', '<f8', ()),
             ('SL12I', '<f8', ()), ('SL22R', '<f8', ()), ('SL22I', '<f8', ()), ('FP2R', '<f8', ()), ('FP2I', '<f8', ()),
             ('FB2R', '<f8', ()), ('FB2I', '<f8', ()), ('EX2R', '<f8', ()), ('EX2I', '<f8', ()), ('EY2R', '<f8', ()),
             ('EY2I', '<f8', ()), ('ET2R', '<f8', ()), ('ET2I', '<f8', ()), ('EL12R', '<f8', ()), ('EL12I', '<f8', ()),
             ('EL22R', '<f8', ()), ('EL22I', '<f8', ()), ('SRP2R', '<f8', ()), ('SRP2I', '<f8', ()),
             ('SRB2R', '<f8', ()), ('SRB2I', '<f8', ()), ('ELTYPE', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FAST(object):
    name = 'FAST'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TX', '<f8', ()), ('TY', '<f8', ()), ('TZ', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FAST_CPLX(object):
    name = 'FAST_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TXR', '<f8', ()), ('TYR', '<f8', ()), ('TZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('TXI', '<f8', ()), ('TYI', '<f8', ()), ('TZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GAP(object):
    name = 'GAP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FX', '<f8', ()), ('SFY', '<f8', ()), ('SFZ', '<f8', ()), ('U', '<f8', ()),
             ('V', '<f8', ()), ('W', '<f8', ()), ('SV', '<f8', ()), ('SW', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GAP_NL(object):
    name = 'GAP_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CPX', '<f8', ()), ('SHY', '<f8', ()), ('SHZ', '<f8', ()), ('AU', '<f8', ()),
             ('SHV', '<f8', ()), ('SHW', '<f8', ()), ('SLV', '<f8', ()), ('SLP', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA(object):
    name = 'HEXA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('X', '<f8', ()),
                                                                                                  ('Y', '<f8', ()),
                                                                                                  ('Z', '<f8', ()),
                                                                                                  ('TXY', '<f8', ()),
                                                                                                  ('TYZ', '<f8', ()),
                                                                                                  ('TZX', '<f8', ())],
                                                                                                 (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA20_27FDNL(object):
    name = 'HEXA20_27FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (27,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA20_8FDNL(object):
    name = 'HEXA20_8FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (8,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA20_FD(object):
    name = 'HEXA20_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      (27,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA_CPLX(object):
    name = 'HEXA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('XR', '<f8', ()),
                                                                                                  ('YR', '<f8', ()),
                                                                                                  ('ZR', '<f8', ()),
                                                                                                  ('TXYR', '<f8', ()),
                                                                                                  ('TYZR', '<f8', ()),
                                                                                                  ('TZXR', '<f8', ()),
                                                                                                  ('XI', '<f8', ()),
                                                                                                  ('YI', '<f8', ()),
                                                                                                  ('ZI', '<f8', ()),
                                                                                                  ('TXYI', '<f8', ()),
                                                                                                  ('TYZI', '<f8', ()),
                                                                                                  ('TZXI', '<f8', ())],
                                                                                                 (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA_FD(object):
    name = 'HEXA_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA_FDNL(object):
    name = 'HEXA_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (8,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA_NL(object):
    name = 'HEXA_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('SX', '<f8', ()),
                                                                                                  ('SY', '<f8', ()),
                                                                                                  ('SZ', '<f8', ()),
                                                                                                  ('SXY', '<f8', ()),
                                                                                                  ('SYZ', '<f8', ()),
                                                                                                  ('SZX', '<f8', ()),
                                                                                                  ('SE', '<f8', ()),
                                                                                                  ('EPS', '<f8', ()),
                                                                                                  ('ECS', '<f8', ()),
                                                                                                  ('EX', '<f8', ()),
                                                                                                  ('EY', '<f8', ()),
                                                                                                  ('EZ', '<f8', ()),
                                                                                                  ('EXY', '<f8', ()),
                                                                                                  ('EYZ', '<f8', ()),
                                                                                                  ('EZX', '<f8', ())],
                                                                                                 (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFHEXA(object):
    name = 'IFHEXA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('GRID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('SXZ', '<f8', ()), ('NX', '<f8', ()), ('NXY', '<f8', ()),
                                                             ('NXZ', '<f8', ()), ('DAMAGE', '<f8', ())], (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFPENTA(object):
    name = 'IFPENTA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('GRID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('SXZ', '<f8', ()), ('NX', '<f8', ()), ('NXY', '<f8', ()),
                                                             ('NXZ', '<f8', ()), ('DAMAGE', '<f8', ())], (7,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA(object):
    name = 'PENTA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('X', '<f8', ()),
                                                                                                  ('Y', '<f8', ()),
                                                                                                  ('Z', '<f8', ()),
                                                                                                  ('TXY', '<f8', ()),
                                                                                                  ('TYZ', '<f8', ()),
                                                                                                  ('TZX', '<f8', ())],
                                                                                                 (7,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA15_21FDNL(object):
    name = 'PENTA15_21FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (21,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA15_6FDNL(object):
    name = 'PENTA15_6FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (6,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA15_FD(object):
    name = 'PENTA15_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      (21,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA_CPLX(object):
    name = 'PENTA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('XR', '<f8', ()),
                                                                                                  ('YR', '<f8', ()),
                                                                                                  ('ZR', '<f8', ()),
                                                                                                  ('TXYR', '<f8', ()),
                                                                                                  ('TYZR', '<f8', ()),
                                                                                                  ('TZXR', '<f8', ()),
                                                                                                  ('XI', '<f8', ()),
                                                                                                  ('YI', '<f8', ()),
                                                                                                  ('ZI', '<f8', ()),
                                                                                                  ('TXYI', '<f8', ()),
                                                                                                  ('TYZI', '<f8', ()),
                                                                                                  ('TZXI', '<f8', ())],
                                                                                                 (7,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA_FD(object):
    name = 'PENTA_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      (6,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA_FDNL(object):
    name = 'PENTA_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (6,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA_NL(object):
    name = 'PENTA_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('SX', '<f8', ()),
                                                                                                  ('SY', '<f8', ()),
                                                                                                  ('SZ', '<f8', ()),
                                                                                                  ('SXY', '<f8', ()),
                                                                                                  ('SYZ', '<f8', ()),
                                                                                                  ('SZX', '<f8', ()),
                                                                                                  ('SE', '<f8', ()),
                                                                                                  ('EPS', '<f8', ()),
                                                                                                  ('ECS', '<f8', ()),
                                                                                                  ('EX', '<f8', ()),
                                                                                                  ('EY', '<f8', ()),
                                                                                                  ('EZ', '<f8', ()),
                                                                                                  ('EXY', '<f8', ()),
                                                                                                  ('EYZ', '<f8', ()),
                                                                                                  ('EZX', '<f8', ())],
                                                                                                 (7,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4(object):
    name = 'QUAD4'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FD1', '<f8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('XY1', '<f8', ()),
             ('FD2', '<f8', ()), ('X2', '<f8', ()), ('Y2', '<f8', ()), ('XY2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_COMP(object):
    name = 'QUAD4_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('PLY', '<i8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('T1', '<f8', ()),
             ('L1', '<f8', ()), ('L2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_COMP_CPLX(object):
    name = 'QUAD4_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('PLY', '<i8', ()), ('X1R', '<f8', ()), ('Y1R', '<f8', ()), ('T1R', '<f8', ()),
             ('L1R', '<f8', ()), ('L2R', '<f8', ()), ('X1I', '<f8', ()), ('Y1I', '<f8', ()), ('T1I', '<f8', ()),
             ('L1I', '<f8', ()), ('L2I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_CPLX(object):
    name = 'QUAD4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FD1', '<f8', ()), ('X1R', '<f8', ()), ('X1I', '<f8', ()), ('Y1R', '<f8', ()),
             ('Y1I', '<f8', ()), ('TXY1R', '<f8', ()), ('TXY1I', '<f8', ()), ('FD2', '<f8', ()), ('X2R', '<f8', ()),
             ('X2I', '<f8', ()), ('Y2R', '<f8', ()), ('Y2I', '<f8', ()), ('TXY2R', '<f8', ()), ('TXY2I', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_FD(object):
    name = 'QUAD4_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()),
             ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()), ('SXY', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_FDNL(object):
    name = 'QUAD4_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PRESSURE', '<f8', ()), ('VOLSTR', '<f8', ()),
                                                             ('EX', '<f8', ()), ('EY', '<f8', ()), ('EZ', '<f8', ()),
                                                             ('EXY', '<f8', ())], (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_FD_CPLX(object):
    name = 'QUAD4_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SXR', '<f8', ()), ('SXI', '<f8', ()),
                                                             ('SYR', '<f8', ()), ('SYI', '<f8', ()),
                                                             ('SXYR', '<f8', ()), ('SXYI', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_NL(object):
    name = 'QUAD4_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('FD1', '<f8', ()), ('SX1', '<f8', ()), ('SY1', '<f8', ()), ('TXY1', '<f8', ()),
             ('ES1', '<f8', ()), ('EPS1', '<f8', ()), ('ECS1', '<f8', ()), ('EX1', '<f8', ()), ('EY1', '<f8', ()),
             ('ETXY1', '<f8', ()), ('FD2', '<f8', ()), ('SX2', '<f8', ()), ('SY2', '<f8', ()), ('TXY2', '<f8', ()),
             ('ES2', '<f8', ()), ('EPS2', '<f8', ()), ('ECS2', '<f8', ()), ('EX2', '<f8', ()), ('EY2', '<f8', ()),
             ('ETXY2', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8(object):
    name = 'QUAD8'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN'
    subtables = []


@register_table
class QUAD8_4FDNL(object):
    name = 'QUAD8_4FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FDNL'
    subtables = []


@register_table
class QUAD8_9FDNL(object):
    name = 'QUAD8_9FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PRESSURE', '<f8', ()), ('VOLSTR', '<f8', ()),
                                                             ('EX', '<f8', ()), ('EY', '<f8', ()), ('EZ', '<f8', ()),
                                                             ('EXY', '<f8', ())], (9,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8_COMP(object):
    name = 'QUAD8_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class QUAD8_COMP_CPLX(object):
    name = 'QUAD8_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class QUAD8_CPLX(object):
    name = 'QUAD8_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN_CPLX'
    subtables = []


@register_table
class QUAD8_FD(object):
    name = 'QUAD8_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()),
             ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()), ('SXY', '<f8', ())], (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8_FD_CPLX(object):
    name = 'QUAD8_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SXR', '<f8', ()), ('SXI', '<f8', ()),
                                                             ('SYR', '<f8', ()), ('SYI', '<f8', ()),
                                                             ('SXYR', '<f8', ()), ('SXYI', '<f8', ())], (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUADR(object):
    name = 'QUADR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN'
    subtables = []


@register_table
class QUADR_COMP(object):
    name = 'QUADR_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class QUADR_COMP_CPLX(object):
    name = 'QUADR_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class QUADR_CPLX(object):
    name = 'QUADR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD_CN_CPLX'
    subtables = []


@register_table
class QUADR_FD(object):
    name = 'QUADR_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3'
    subtables = []


@register_table
class QUADR_FD_CPLX(object):
    name = 'QUADR_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_CPLX'
    subtables = []


@register_table
class QUADR_NL(object):
    name = 'QUADR_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_NL'
    subtables = []


@register_table
class QUADX4_FD(object):
    name = 'QUADX4_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FD'
    subtables = []


@register_table
class QUADX4_FDNL(object):
    name = 'QUADX4_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FDNL'
    subtables = []


@register_table
class QUADX4_FD_CPLX(object):
    name = 'QUADX4_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FD_CPLX'
    subtables = []


@register_table
class QUADX8_4FDNL(object):
    name = 'QUADX8_4FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_FDNL'
    subtables = []


@register_table
class QUADX8_9FDNL(object):
    name = 'QUADX8_9FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_9FDNL'
    subtables = []


@register_table
class QUADX8_FD(object):
    name = 'QUADX8_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_FD'
    subtables = []


@register_table
class QUADX8_FD_CPLX(object):
    name = 'QUADX8_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD8_FD_CPLX'
    subtables = []


@register_table
class QUAD_CN(object):
    name = 'QUAD_CN'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('SS', [('GRID', '<i8', ()), ('FD1', '<f8', ()), ('X1', '<f8', ()),
                                                             ('Y1', '<f8', ()), ('TXY1', '<f8', ()), ('FD2', '<f8', ()),
                                                             ('X2', '<f8', ()), ('Y2', '<f8', ()), ('TXY2', '<f8', ())],
                                                      (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD_CN_CPLX(object):
    name = 'QUAD_CN_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('SS',
                                                      [('GRID', '<i8', ()), ('FD1', '<f8', ()), ('X1R', '<f8', ()),
                                                       ('X1I', '<f8', ()), ('Y1R', '<f8', ()), ('Y1I', '<f8', ()),
                                                       ('TXY1R', '<f8', ()), ('TXY1I', '<f8', ()), ('FD2', '<f8', ()),
                                                       ('X2R', '<f8', ()), ('X2I', '<f8', ()), ('Y2R', '<f8', ()),
                                                       ('Y2I', '<f8', ()), ('TXY2R', '<f8', ()), ('TXY2I', '<f8', ())],
                                                      (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RAC2D(object):
    name = 'RAC2D'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
             ('TXY', '<f8', ()), ('KI', '<f8', ()), ('KII', '<f8', ()), ('S8', '<f8', ()), ('S9', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RAC3D(object):
    name = 'RAC3D'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('TXY', '<f8', ()),
             ('TYZ', '<f8', ()), ('TZX', '<f8', ()), ('KI', '<f8', ()), ('KII', '<f8', ()), ('KIII', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROD(object):
    name = 'ROD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('A', '<f8', ()), ('MSA', '<f8', ()), ('T', '<f8', ()), ('MST', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROD_CPLX(object):
    name = 'ROD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AR', '<f8', ()), ('AI', '<f8', ()), ('TR', '<f8', ()), ('TI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROD_NL(object):
    name = 'ROD_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('SE', '<f8', ()), ('TE', '<f8', ()), ('EPS', '<f8', ()),
             ('ECS', '<f8', ()), ('LTS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ROD_RR(object):
    name = 'ROD_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('A', '<f8', ()), ('T', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SEAM(object):
    name = 'SEAM'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('BMAXA', '<f8', ()), ('BMINA', '<f8', ()), ('BMAXB', '<f8', ()),
             ('BMINB', '<f8', ()), ('TMAX', '<f8', ()), ('BRNG', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SEAMP(object):
    name = 'SEAMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/HEXA'
    subtables = []


@register_table
class SEAMP_CPLX(object):
    name = 'SEAMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/HEXA_CPLX'
    subtables = []


@register_table
class SHEAR(object):
    name = 'SHEAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TMAX', '<f8', ()), ('TAVG', '<f8', ()), ('MS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SHEAR_CPLX(object):
    name = 'SHEAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TMAXR', '<f8', ()), ('TMAXI', '<f8', ()), ('TAVGR', '<f8', ()), ('TAVGI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SHEAR_RR(object):
    name = 'SHEAR_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TMAX', '<f8', ()), ('TAVG', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOT3(object):
    name = 'SLOT3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RC', '<f8', ()), ('AC', '<f8', ()), ('TE1', '<f8', ()), ('TE2', '<f8', ()),
             ('TE3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOT3_CPLX(object):
    name = 'SLOT3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RCR', '<f8', ()), ('ACR', '<f8', ()), ('TE1R', '<f8', ()), ('TE2R', '<f8', ()),
             ('TE3R', '<f8', ()), ('RCI', '<f8', ()), ('ACI', '<f8', ()), ('TE1I', '<f8', ()), ('TE2I', '<f8', ()),
             ('TE3I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOT4(object):
    name = 'SLOT4'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RC', '<f8', ()), ('AC', '<f8', ()), ('TE1', '<f8', ()), ('TE2', '<f8', ()),
             ('TE3', '<f8', ()), ('TE4', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SLOT4_CPLX(object):
    name = 'SLOT4_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('RCR', '<f8', ()), ('ACR', '<f8', ()), ('TE1R', '<f8', ()), ('TE2R', '<f8', ()),
             ('TE3R', '<f8', ()), ('TE4R', '<f8', ()), ('RCI', '<f8', ()), ('ACI', '<f8', ()), ('TE1I', '<f8', ()),
             ('TE2I', '<f8', ()), ('TE3I', '<f8', ()), ('TE4I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA(object):
    name = 'TETRA'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('X', '<f8', ()),
                                                                                                  ('Y', '<f8', ()),
                                                                                                  ('Z', '<f8', ()),
                                                                                                  ('TXY', '<f8', ()),
                                                                                                  ('TYZ', '<f8', ()),
                                                                                                  ('TZX', '<f8', ())],
                                                                                                 (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA10_4FDNL(object):
    name = 'TETRA10_4FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA10_5FDNL(object):
    name = 'TETRA10_5FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA10_FD(object):
    name = 'TETRA10_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      (5,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA4_FDNL(object):
    name = 'TETRA4_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA_CPLX(object):
    name = 'TETRA_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('XR', '<f8', ()),
                                                                                                  ('YR', '<f8', ()),
                                                                                                  ('ZR', '<f8', ()),
                                                                                                  ('TXYR', '<f8', ()),
                                                                                                  ('TYZR', '<f8', ()),
                                                                                                  ('TZXR', '<f8', ()),
                                                                                                  ('XI', '<f8', ()),
                                                                                                  ('YI', '<f8', ()),
                                                                                                  ('ZI', '<f8', ()),
                                                                                                  ('TXYI', '<f8', ()),
                                                                                                  ('TYZI', '<f8', ()),
                                                                                                  ('TZXI', '<f8', ())],
                                                                                                 (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA_FD(object):
    name = 'TETRA_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PA', '<f8', ()), ('AX', '<f8', ()), ('AY', '<f8', ()),
                                                             ('AZ', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('SY', '<f8', ()), ('SYZ', '<f8', ()), ('PB', '<f8', ()),
                                                             ('BX', '<f8', ()), ('BY', '<f8', ()), ('BZ', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SZX', '<f8', ()), ('PC', '<f8', ()),
                                                             ('CX', '<f8', ()), ('CY', '<f8', ()), ('CZ', '<f8', ())],
                                                      ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA_FDNL(object):
    name = 'TETRA_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()), ('SYZ', '<f8', ()),
                                                             ('SZX', '<f8', ()), ('PRESSURE', '<f8', ()),
                                                             ('VOLSTR', '<f8', ()), ('EX', '<f8', ()),
                                                             ('EY', '<f8', ()), ('EZ', '<f8', ()), ('EXY', '<f8', ()),
                                                             ('EYZ', '<f8', ()), ('EZX', '<f8', ())], ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TETRA_NL(object):
    name = 'TETRA_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('CID', '<i8', ()), ('CTYPE', 'S4', ()), ('NODEF', '<i8', ()), ('SS',
                                                                                                 [('GRID', '<i8', ()),
                                                                                                  ('SX', '<f8', ()),
                                                                                                  ('SY', '<f8', ()),
                                                                                                  ('SZ', '<f8', ()),
                                                                                                  ('SXY', '<f8', ()),
                                                                                                  ('SYZ', '<f8', ()),
                                                                                                  ('SZX', '<f8', ()),
                                                                                                  ('SE', '<f8', ()),
                                                                                                  ('EPS', '<f8', ()),
                                                                                                  ('ECS', '<f8', ()),
                                                                                                  ('EX', '<f8', ()),
                                                                                                  ('EY', '<f8', ()),
                                                                                                  ('EZ', '<f8', ()),
                                                                                                  ('EXY', '<f8', ()),
                                                                                                  ('EYZ', '<f8', ()),
                                                                                                  ('EZX', '<f8', ())],
                                                                                                 (5,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3(object):
    name = 'TRIA3'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('FD1', '<f8', ()), ('X1', '<f8', ()), ('Y1', '<f8', ()), ('TXY1', '<f8', ()),
                                         ('FD2', '<f8', ()), ('X2', '<f8', ()), ('Y2', '<f8', ()), ('TXY2', '<f8', ())],
                                  ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_1FDNL(object):
    name = 'TRIA3_1FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PRESSURE', '<f8', ()), ('VOLSTR', '<f8', ()),
                                                             ('EX', '<f8', ()), ('EY', '<f8', ()), ('EZ', '<f8', ()),
                                                             ('EXY', '<f8', ())], ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_3FDNL(object):
    name = 'TRIA3_3FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()),
                                                             ('SZ', '<f8', ()), ('SXY', '<f8', ()),
                                                             ('PRESSURE', '<f8', ()), ('VOLSTR', '<f8', ()),
                                                             ('EX', '<f8', ()), ('EY', '<f8', ()), ('EZ', '<f8', ()),
                                                             ('EXY', '<f8', ())], (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_COMP(object):
    name = 'TRIA3_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class TRIA3_COMP_CPLX(object):
    name = 'TRIA3_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class TRIA3_CPLX(object):
    name = 'TRIA3_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS', [('FD1', '<f8', ()), ('X1R', '<f8', ()), ('X1I', '<f8', ()), ('Y1R', '<f8', ()),
                                         ('Y1I', '<f8', ()), ('TXY1R', '<f8', ()), ('TXY1I', '<f8', ()),
                                         ('FD2', '<f8', ()), ('X2R', '<f8', ()), ('X2I', '<f8', ()), ('Y2R', '<f8', ()),
                                         ('Y2I', '<f8', ()), ('TXY2R', '<f8', ()), ('TXY2I', '<f8', ())], ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_FD(object):
    name = 'TRIA3_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()),
             ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()), ('SXY', '<f8', ())], ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_FD_CPLX(object):
    name = 'TRIA3_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SXR', '<f8', ()), ('SXI', '<f8', ()),
                                                             ('SYR', '<f8', ()), ('SYI', '<f8', ()),
                                                             ('SXYR', '<f8', ()), ('SXYI', '<f8', ())], ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_NL(object):
    name = 'TRIA3_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_NL'
    subtables = []


@register_table
class TRIA6(object):
    name = 'TRIA6'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('SS', [('GRID', '<i8', ()), ('FD1', '<f8', ()), ('X1', '<f8', ()),
                                                             ('Y1', '<f8', ()), ('TXY1', '<f8', ()), ('FD2', '<f8', ()),
                                                             ('X2', '<f8', ()), ('Y2', '<f8', ()), ('TXY2', '<f8', ())],
                                                      (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA6_COMP(object):
    name = 'TRIA6_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class TRIA6_COMP_CPLX(object):
    name = 'TRIA6_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class TRIA6_CPLX(object):
    name = 'TRIA6_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TERM', 'S4', ()), ('SS',
                                                      [('GRID', '<i8', ()), ('FD1', '<f8', ()), ('X1R', '<f8', ()),
                                                       ('X1I', '<f8', ()), ('Y1R', '<f8', ()), ('Y1I', '<f8', ()),
                                                       ('TXY1R', '<f8', ()), ('TXY1I', '<f8', ()), ('FD2', '<f8', ()),
                                                       ('X2R', '<f8', ()), ('X2I', '<f8', ()), ('Y2R', '<f8', ()),
                                                       ('Y2I', '<f8', ()), ('TXY2R', '<f8', ()), ('TXY2I', '<f8', ())],
                                                      (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA6_FD(object):
    name = 'TRIA6_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()),
             ('SS', [('ID', '<i8', ()), ('SX', '<f8', ()), ('SY', '<f8', ()), ('SXY', '<f8', ())], (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA6_FDNL(object):
    name = 'TRIA6_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_3FDNL'
    subtables = []


@register_table
class TRIA6_FD_CPLX(object):
    name = 'TRIA6_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('TYPE', 'S4', ()), ('SS', [('ID', '<i8', ()), ('SXR', '<f8', ()), ('SXI', '<f8', ()),
                                                             ('SYR', '<f8', ()), ('SYI', '<f8', ()),
                                                             ('SXYR', '<f8', ()), ('SXYI', '<f8', ())], (3,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIAR(object):
    name = 'TRIAR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6'
    subtables = []


@register_table
class TRIAR_1FD(object):
    name = 'TRIAR_1FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3'
    subtables = []


@register_table
class TRIAR_1FD_CPLX(object):
    name = 'TRIAR_1FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_CPLX'
    subtables = []


@register_table
class TRIAR_4FD(object):
    name = 'TRIAR_4FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6'
    subtables = []


@register_table
class TRIAR_4FD_CPLX(object):
    name = 'TRIAR_4FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_CPLX'
    subtables = []


@register_table
class TRIAR_COMP(object):
    name = 'TRIAR_COMP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP'
    subtables = []


@register_table
class TRIAR_COMP_CPLX(object):
    name = 'TRIAR_COMP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX'
    subtables = []


@register_table
class TRIAR_CPLX(object):
    name = 'TRIAR_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_CPLX'
    subtables = []


@register_table
class TRIAR_NL(object):
    name = 'TRIAR_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/QUAD4_NL'
    subtables = []


@register_table
class TRIAX3_1FDNL(object):
    name = 'TRIAX3_1FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_1FDNL'
    subtables = []


@register_table
class TRIAX3_3FDNL(object):
    name = 'TRIAX3_3FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_3FDNL'
    subtables = []


@register_table
class TRIAX3_FD(object):
    name = 'TRIAX3_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_FD'
    subtables = []


@register_table
class TRIAX3_FD_CPLX(object):
    name = 'TRIAX3_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_FD_CPLX'
    subtables = []


@register_table
class TRIAX6(object):
    name = 'TRIAX6'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), (
    'SS', [('LOC', '<i8', ()), ('RS', '<f8', ()), ('AZS', '<f8', ()), ('AS', '<f8', ()), ('TS', '<f8', ())], (4,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIAX6_CPLX(object):
    name = 'TRIAX6_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('LOC', '<i8', ()), ('RSR', '<f8', ()), ('RSI', '<f8', ()), ('AZSR', '<f8', ()),
                                   ('AZSI', '<f8', ()), ('ASR', '<f8', ()), ('ASI', '<f8', ()), ('TSR', '<f8', ()),
                                   ('TSI', '<f8', ())], (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIAX6_FD(object):
    name = 'TRIAX6_FD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_FD'
    subtables = []


@register_table
class TRIAX6_FDNL(object):
    name = 'TRIAX6_FDNL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA3_3FDNL'
    subtables = []


@register_table
class TRIAX6_FD_CPLX(object):
    name = 'TRIAX6_FD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/TRIA6_FD_CPLX'
    subtables = []


@register_table
class TUBE(object):
    name = 'TUBE'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('MSA', '<f8', ()), ('TS', '<f8', ()), ('MST', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TUBE_CPLX(object):
    name = 'TUBE_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('ASR', '<f8', ()), ('ASI', '<f8', ()), ('TSR', '<f8', ()), ('TSI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TUBE_NL(object):
    name = 'TUBE_NL'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('SE', '<f8', ()), ('TE', '<f8', ()), ('EPS', '<f8', ()),
             ('ECS', '<f8', ()), ('LTS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TUBE_RR(object):
    name = 'TUBE_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('TS', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VISC_CPLX(object):
    name = 'VISC_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('ASR', '<f8', ()), ('ASI', '<f8', ()), ('TAUR', '<f8', ()), ('TAUI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VISC_RR(object):
    name = 'VISC_RR'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('TAU', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WELD(object):
    name = 'WELD'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('AS', '<f8', ()), ('MXA', '<f8', ()), ('MNA', '<f8', ()), ('MXB', '<f8', ()),
             ('MNB', '<f8', ()), ('TMX', '<f8', ()), ('BRNG', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WELD_CPLX(object):
    name = 'WELD_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = [('EID', '<i8', ()), ('ASR', '<f8', ()), ('ASI', '<f8', ()), ('MXAR', '<f8', ()), ('MXAI', '<f8', ()),
             ('MNAR', '<f8', ()), ('MNAI', '<f8', ()), ('MXBR', '<f8', ()), ('MXBI', '<f8', ()), ('MNBR', '<f8', ()),
             ('MNBI', '<f8', ()), ('TMXR', '<f8', ()), ('TMXI', '<f8', ()), ('BRNGR', '<f8', ()), ('BRNGI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WELDC(object):
    name = 'WELDC'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD'
    subtables = []


@register_table
class WELDC_CPLX(object):
    name = 'WELDC_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD_CPLX'
    subtables = []


@register_table
class WELDP(object):
    name = 'WELDP'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD'
    subtables = []


@register_table
class WELDP_CPLX(object):
    name = 'WELDP_CPLX'
    path = '/NASTRAN/RESULT/ELEMENTAL/STRESS'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/ELEMENTAL/STRESS/WELD_CPLX'
    subtables = []


@register_table
class VIBRATION_INTENSITY(object):
    name = 'VIBRATION_INTENSITY'
    path = '/NASTRAN/RESULT/ELEMENTAL'
    dtype = [('EID', '<i8', ()), ('ELTYPE', '<i8', ()), ('UNX', '<f8', ()), ('UNY', '<f8', ()), ('UNZ', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPOTWFTG(object):
    name = 'SPOTWFTG'
    path = '/NASTRAN/RESULT/FATIGUE/CBAR'
    dtype = [('GRID', '<i8', ()), ('ANGL', '<f8', ()), ('LIFE', '<f8', ()), ('LOGLIFE', '<f8', ()),
             ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()), ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()),
             ('maxforce', '<f8', ()), ('MAXSTRES', '<f8', ()), ('MINSTRES', '<f8', ()), ('FOS', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/FATIGUE/CBAR'
    dtype = [('EID', '<i8', ()), ('SPOTWFTG_POS', '<i8', ()), ('SPOTWFTG_LEN', '<i8', ()), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/FATIGUE/CBAR/SPOTWFTG']


@register_table
class CHEXA(object):
    name = 'CHEXA'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('GRID', '<i8', ()), ('LAYER', '<i8', ()), ('LIFE', '<f8', ()), ('LOGLIFE', '<f8', ()),
             ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()), ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()),
             ('MAXVALUE', '<f8', ()), ('MINVALUE', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXA_CN(object):
    name = 'CHEXA_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER', '<i8', ()), ('LIFE', '<f8', ()),
                                         ('LOGLIFE', '<f8', ()), ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()),
                                         ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()), ('MAXVALUE', '<f8', ()),
                                         ('MINVALUE', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ())], (8,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CPENTA(object):
    name = 'CPENTA'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CHEXA'
    subtables = []


@register_table
class CPENTA_CN(object):
    name = 'CPENTA_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER', '<i8', ()), ('LIFE', '<f8', ()),
                                         ('LOGLIFE', '<f8', ()), ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()),
                                         ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()), ('MAXVALUE', '<f8', ()),
                                         ('MINVALUE', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ())], (6,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4(object):
    name = 'CQUAD4'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('LIFE1', '<f8', ()),
             ('LOGLIFE1', '<f8', ()), ('EQLIFE1', '<f8', ()), ('LGEQLIF1', '<f8', ()), ('DAMAGE1', '<f8', ()),
             ('LGDAMAG1', '<f8', ()), ('MAXVALUE1', '<f8', ()), ('MINVALUE1', '<f8', ()), ('FOS1', '<f8', ()),
             ('CRTANGL1', '<f8', ()), ('LAYER2', '<i8', ()), ('LIFE2', '<f8', ()), ('LOGLIFE2', '<f8', ()),
             ('EQLIFE2', '<f8', ()), ('LGEQLIF2', '<f8', ()), ('DAMAGE2', '<f8', ()), ('LGDAMAG2', '<f8', ()),
             ('MAXVALUE2', '<f8', ()), ('MINVALUE2', '<f8', ()), ('FOS2', '<f8', ()), ('CRTANGL2', '<f8', ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_AUTO(object):
    name = 'CQUAD4_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()),
                                         ('PLANG2', '<f8', ()), ('DSTS2', '<f8', ())], ()), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_AUTO_CN(object):
    name = 'CQUAD4_BM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()),
                                         ('PLANG2', '<f8', ()), ('DSTS2', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_SIMPLE(object):
    name = 'CQUAD4_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('BIAXS2', '<f8', ()), ('ANGLS2', '<f8', ()),
                                         ('ANGLM2', '<f8', ()), ('BIAXG2', '<f8', ()), ('MAXAP2', '<f8', ())], ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_SIMPLE_CN(object):
    name = 'CQUAD4_BM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('BIAXS2', '<f8', ()), ('ANGLS2', '<f8', ()),
                                         ('ANGLM2', '<f8', ()), ('BIAXG2', '<f8', ()), ('MAXAP2', '<f8', ())], (4,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_SM_AUTO_CN(object):
    name = 'CQUAD4_BM_SM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ())], (4,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_SM_SIMPLE_CN(object):
    name = 'CQUAD4_BM_SM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_SM_STD_CN(object):
    name = 'CQUAD4_BM_SM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_STD(object):
    name = 'CQUAD4_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()), ('DSTS2', '<f8', ())], ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_BM_STD_CN(object):
    name = 'CQUAD4_BM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()), ('DSTS2', '<f8', ())], (4,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_CN(object):
    name = 'CQUAD4_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('LIFE1', '<f8', ()),
                                         ('LOGLIFE1', '<f8', ()), ('LGEQLIF1', '<f8', ()), ('EQLIFE1', '<f8', ()),
                                         ('DAMAGE1', '<f8', ()), ('LGDAMAG1', '<f8', ()), ('MAXVALUE1', '<f8', ()),
                                         ('MINVALUE1', '<f8', ()), ('FOS1', '<f8', ()), ('CRTANGL1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('LIFE2', '<f8', ()), ('LOGLIFE2', '<f8', ()),
                                         ('EQLIFE2', '<f8', ()), ('LGEQLIF2', '<f8', ()), ('DAMAGE2', '<f8', ()),
                                         ('LGDAMAG2', '<f8', ()), ('MAXVALUE2', '<f8', ()), ('MINVALUE2', '<f8', ()),
                                         ('FOS2', '<f8', ()), ('CRTANGL2', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_SM(object):
    name = 'CQUAD4_SM'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('LIFE', '<f8', ()),
                                         ('LOGLIFE', '<f8', ()), ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()),
                                         ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()), ('MAXSTRES', '<f8', ()),
                                         ('MINSTRES', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ())], (4,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD8(object):
    name = 'CQUAD8'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4'
    subtables = []


@register_table
class CQUAD8_BM_AUTO(object):
    name = 'CQUAD8_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_AUTO'
    subtables = []


@register_table
class CQUAD8_BM_AUTO_CN(object):
    name = 'CQUAD8_BM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_AUTO_CN'
    subtables = []


@register_table
class CQUAD8_BM_SIMPLE(object):
    name = 'CQUAD8_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SIMPLE'
    subtables = []


@register_table
class CQUAD8_BM_SIMPLE_CN(object):
    name = 'CQUAD8_BM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SIMPLE_CN'
    subtables = []


@register_table
class CQUAD8_BM_SM_AUTO_CN(object):
    name = 'CQUAD8_BM_SM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_AUTO_CN'
    subtables = []


@register_table
class CQUAD8_BM_SM_SIMPLE_CN(object):
    name = 'CQUAD8_BM_SM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_SIMPLE_CN'
    subtables = []


@register_table
class CQUAD8_BM_SM_STD_CN(object):
    name = 'CQUAD8_BM_SM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_STD_CN'
    subtables = []


@register_table
class CQUAD8_BM_STD(object):
    name = 'CQUAD8_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_STD'
    subtables = []


@register_table
class CQUAD8_BM_STD_CN(object):
    name = 'CQUAD8_BM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_STD_CN'
    subtables = []


@register_table
class CQUAD8_CN(object):
    name = 'CQUAD8_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_CN'
    subtables = []


@register_table
class CQUAD8_SM(object):
    name = 'CQUAD8_SM'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_SM'
    subtables = []


@register_table
class CQUADR(object):
    name = 'CQUADR'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4'
    subtables = []


@register_table
class CQUADR_BM_AUTO(object):
    name = 'CQUADR_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_AUTO'
    subtables = []


@register_table
class CQUADR_BM_AUTO_CN(object):
    name = 'CQUADR_BM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_AUTO_CN'
    subtables = []


@register_table
class CQUADR_BM_SIMPLE(object):
    name = 'CQUADR_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SIMPLE'
    subtables = []


@register_table
class CQUADR_BM_SIMPLE_CN(object):
    name = 'CQUADR_BM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SIMPLE_CN'
    subtables = []


@register_table
class CQUADR_BM_SM_AUTO_CN(object):
    name = 'CQUADR_BM_SM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_AUTO_CN'
    subtables = []


@register_table
class CQUADR_BM_SM_SIMPLE_CN(object):
    name = 'CQUADR_BM_SM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_SIMPLE_CN'
    subtables = []


@register_table
class CQUADR_BM_SM_STD_CN(object):
    name = 'CQUADR_BM_SM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_SM_STD_CN'
    subtables = []


@register_table
class CQUADR_BM_STD(object):
    name = 'CQUADR_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_STD'
    subtables = []


@register_table
class CQUADR_BM_STD_CN(object):
    name = 'CQUADR_BM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_BM_STD_CN'
    subtables = []


@register_table
class CQUADR_CN(object):
    name = 'CQUADR_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_CN'
    subtables = []


@register_table
class CQUADR_SM(object):
    name = 'CQUADR_SM'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4_SM'
    subtables = []


@register_table
class CSHEAR(object):
    name = 'CSHEAR'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CHEXA'
    subtables = []


@register_table
class CTETRA(object):
    name = 'CTETRA'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CHEXA'
    subtables = []


@register_table
class CTETRA_CN(object):
    name = 'CTETRA_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER', '<i8', ()), ('LIFE', '<f8', ()),
                                         ('LOGLIFE', '<f8', ()), ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()),
                                         ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()), ('MAXVALUE', '<f8', ()),
                                         ('MINVALUE', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ())], (4,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA3(object):
    name = 'CTRIA3'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4'
    subtables = []


@register_table
class CTRIA3_BM_AUTO(object):
    name = 'CTRIA3_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_AUTO'
    subtables = []


@register_table
class CTRIA3_BM_SIMPLE(object):
    name = 'CTRIA3_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_SIMPLE'
    subtables = []


@register_table
class CTRIA3_BM_STD(object):
    name = 'CTRIA3_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_STD'
    subtables = []


@register_table
class CTRIA6(object):
    name = 'CTRIA6'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4'
    subtables = []


@register_table
class CTRIA6_BM_AUTO(object):
    name = 'CTRIA6_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_AUTO'
    subtables = []


@register_table
class CTRIA6_BM_AUTO_CN(object):
    name = 'CTRIA6_BM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_AUTO_CN'
    subtables = []


@register_table
class CTRIA6_BM_SIMPLE(object):
    name = 'CTRIA6_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_SIMPLE'
    subtables = []


@register_table
class CTRIA6_BM_SIMPLE_CN(object):
    name = 'CTRIA6_BM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_SIMPLE_CN'
    subtables = []


@register_table
class CTRIA6_BM_SM_AUTO_CN(object):
    name = 'CTRIA6_BM_SM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ())], (3,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA6_BM_SM_SIMPLE_CN(object):
    name = 'CTRIA6_BM_SM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ())], (3,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA6_BM_SM_STD_CN(object):
    name = 'CTRIA6_BM_SM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ())], (3,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA6_BM_STD(object):
    name = 'CTRIA6_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_STD'
    subtables = []


@register_table
class CTRIA6_BM_STD_CN(object):
    name = 'CTRIA6_BM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_BM_STD_CN'
    subtables = []


@register_table
class CTRIA6_CN(object):
    name = 'CTRIA6_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_CN'
    subtables = []


@register_table
class CTRIA6_SM(object):
    name = 'CTRIA6_SM'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIAR_SM'
    subtables = []


@register_table
class CTRIAR(object):
    name = 'CTRIAR'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CQUAD4'
    subtables = []


@register_table
class CTRIAR_BM_AUTO(object):
    name = 'CTRIAR_BM_AUTO'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()),
                                         ('PLANG2', '<f8', ()), ('DSTS2', '<f8', ())], ()), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_BM_AUTO_CN(object):
    name = 'CTRIAR_BM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('PLANG1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()),
                                         ('PLANG2', '<f8', ()), ('DSTS2', '<f8', ())], (3,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_BM_SIMPLE(object):
    name = 'CTRIAR_BM_SIMPLE'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('BIAXS2', '<f8', ()), ('ANGLS2', '<f8', ()),
                                         ('ANGLM2', '<f8', ()), ('BIAXG2', '<f8', ()), ('MAXAP2', '<f8', ())], ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_BM_SIMPLE_CN(object):
    name = 'CTRIAR_BM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('BIAXS1', '<f8', ()), ('ANGLS1', '<f8', ()), ('ANGLM1', '<f8', ()),
                                         ('BIAXG1', '<f8', ()), ('MAXAP1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('BIAXS2', '<f8', ()), ('ANGLS2', '<f8', ()),
                                         ('ANGLM2', '<f8', ()), ('BIAXG2', '<f8', ()), ('MAXAP2', '<f8', ())], (3,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_BM_SM_AUTO_CN(object):
    name = 'CTRIAR_BM_SM_AUTO_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIA6_BM_SM_AUTO_CN'
    subtables = []


@register_table
class CTRIAR_BM_SM_SIMPLE_CN(object):
    name = 'CTRIAR_BM_SM_SIMPLE_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIA6_BM_SM_SIMPLE_CN'
    subtables = []


@register_table
class CTRIAR_BM_SM_STD_CN(object):
    name = 'CTRIAR_BM_SM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE/CTRIA6_BM_SM_STD_CN'
    subtables = []


@register_table
class CTRIAR_BM_STD(object):
    name = 'CTRIAR_BM_STD'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()), ('DSTS2', '<f8', ())], ()),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_BM_STD_CN(object):
    name = 'CTRIAR_BM_STD_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('BIAXM1', '<f8', ()),
                                         ('NONPR1', '<f8', ()), ('DSTS1', '<f8', ()), ('LAYER2', '<i8', ()),
                                         ('BIAXM2', '<f8', ()), ('NONPR2', '<f8', ()), ('DSTS2', '<f8', ())], (3,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_CN(object):
    name = 'CTRIAR_CN'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('LIFE1', '<f8', ()),
                                         ('LOGLIFE1', '<f8', ()), ('LGEQLIF1', '<f8', ()), ('EQLIFE1', '<f8', ()),
                                         ('DAMAGE1', '<f8', ()), ('LGDAMAG1', '<f8', ()), ('MAXVALUE1', '<f8', ()),
                                         ('MINVALUE1', '<f8', ()), ('FOS1', '<f8', ()), ('CRTANGL1', '<f8', ()),
                                         ('LAYER2', '<i8', ()), ('LIFE2', '<f8', ()), ('LOGLIFE2', '<f8', ()),
                                         ('EQLIFE2', '<f8', ()), ('LGEQLIF2', '<f8', ()), ('DAMAGE2', '<f8', ()),
                                         ('LGDAMAG2', '<f8', ()), ('MAXVALUE2', '<f8', ()), ('MINVALUE2', '<f8', ()),
                                         ('FOS2', '<f8', ()), ('CRTANGL2', '<f8', ())], (3,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR_SM(object):
    name = 'CTRIAR_SM'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('EID', '<i8', ()), ('SS', [('GRID', '<i8', ()), ('slcty', '<i8', ()), ('LIFE', '<f8', ()),
                                         ('LOGLIFE', '<f8', ()), ('EQLIFE', '<f8', ()), ('LOGEQLIF', '<f8', ()),
                                         ('DAMAGE', '<f8', ()), ('LOGDAMAG', '<f8', ()), ('MAXSTRES', '<f8', ()),
                                         ('MINSTRES', '<f8', ()), ('FOS', '<f8', ()), ('CRITANGL', '<f8', ())], (3,)),
             ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/FATIGUE'
    dtype = [('IDENT_ID', '<i8', ()), ('NANG3', '<i8', ()), ('EVNTNAM', 'S56', ()), ('EQNAME', 'S20', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXA(object):
    name = 'CHEXA'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M0', '<f8', ()), ('M1', '<f8', ()),
                                   ('M2', '<f8', ()), ('M4', '<f8', ()), ('E0', '<f8', ()), ('EP', '<f8', ()),
                                   ('GAMMA', '<f8', ()), ('RMSTRS', '<f8', ()), ('RMSTRN', '<f8', ()),
                                   ('MNSTRS', '<f8', ()), ('MNSTRN', '<f8', ()), ('STRPSGM', '<f8', ()),
                                   ('STRMSGM', '<f8', ()), ('STNPSGM', '<f8', ()), ('STNMSGM', '<f8', ()),
                                   ('PLINDEX', '<f8', ())], ()), ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CHEXA_CN(object):
    name = 'CHEXA_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M0', '<f8', ()), ('M1', '<f8', ()),
                                   ('M2', '<f8', ()), ('M4', '<f8', ()), ('E0', '<f8', ()), ('EP', '<f8', ()),
                                   ('GAMMA', '<f8', ()), ('RMSTRS', '<f8', ()), ('RMSTRN', '<f8', ()),
                                   ('MNSTRS', '<f8', ()), ('MNSTRN', '<f8', ()), ('STRPSGM', '<f8', ()),
                                   ('STRMSGM', '<f8', ()), ('STNPSGM', '<f8', ()), ('STNMSGM', '<f8', ()),
                                   ('PLINDEX', '<f8', ())], (8,)), ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CPENTA(object):
    name = 'CPENTA'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CHEXA'
    subtables = []


@register_table
class CPENTA_CN(object):
    name = 'CPENTA_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M0', '<f8', ()), ('M1', '<f8', ()),
                                   ('M2', '<f8', ()), ('M4', '<f8', ()), ('E0', '<f8', ()), ('EP', '<f8', ()),
                                   ('GAMMA', '<f8', ()), ('RMSTRS', '<f8', ()), ('RMSTRN', '<f8', ()),
                                   ('MNSTRS', '<f8', ()), ('MNSTRN', '<f8', ()), ('STRPSGM', '<f8', ()),
                                   ('STRMSGM', '<f8', ()), ('STNPSGM', '<f8', ()), ('STNMSGM', '<f8', ()),
                                   ('PLINDEX', '<f8', ())], (6,)), ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4(object):
    name = 'CQUAD4'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M01', '<f8', ()), ('M11', '<f8', ()),
                                   ('M21', '<f8', ()), ('M41', '<f8', ()), ('E01', '<f8', ()), ('EP1', '<f8', ()),
                                   ('GAMMA1', '<f8', ()), ('RMSTRS1', '<f8', ()), ('RMSTRN1', '<f8', ()),
                                   ('MNSTRS1', '<f8', ()), ('MNSTRN1', '<f8', ()), ('STRPSGM1', '<f8', ()),
                                   ('STRMSGM1', '<f8', ()), ('STNPSGM1', '<f8', ()), ('STNMSGM1', '<f8', ()),
                                   ('PLINDEX1', '<f8', ()), ('LAYER2', '<i8', ()), ('M02', '<f8', ()),
                                   ('M12', '<f8', ()), ('M22', '<f8', ()), ('M42', '<f8', ()), ('E02', '<f8', ()),
                                   ('EP2', '<f8', ()), ('GAMMA2', '<f8', ()), ('RMSTRS2', '<f8', ()),
                                   ('RMSTRN2', '<f8', ()), ('MNSTRS2', '<f8', ()), ('MNSTRN2', '<f8', ()),
                                   ('STRPSGM2', '<f8', ()), ('STRMSGM2', '<f8', ()), ('STNPSGM2', '<f8', ()),
                                   ('STNMSGM2', '<f8', ()), ('PLINDEX2', '<f8', ())], ()), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD4_CN(object):
    name = 'CQUAD4_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M01', '<f8', ()), ('M11', '<f8', ()),
                                   ('M21', '<f8', ()), ('M41', '<f8', ()), ('E01', '<f8', ()), ('EP1', '<f8', ()),
                                   ('GAMMA1', '<f8', ()), ('RMSTRS1', '<f8', ()), ('RMSTRN1', '<f8', ()),
                                   ('MNSTRS1', '<f8', ()), ('MNSTRN1', '<f8', ()), ('STRPSGM1', '<f8', ()),
                                   ('STRMSGM1', '<f8', ()), ('STNPSGM1', '<f8', ()), ('STNMSGM1', '<f8', ()),
                                   ('PLINDEX1', '<f8', ()), ('LAYER2', '<i8', ()), ('M02', '<f8', ()),
                                   ('M12', '<f8', ()), ('M22', '<f8', ()), ('M42', '<f8', ()), ('E02', '<f8', ()),
                                   ('EP2', '<f8', ()), ('GAMMA2', '<f8', ()), ('RMSTRS2', '<f8', ()),
                                   ('RMSTRN2', '<f8', ()), ('MNSTRS2', '<f8', ()), ('MNSTRN2', '<f8', ()),
                                   ('STRPSGM2', '<f8', ()), ('STRMSGM2', '<f8', ()), ('STNPSGM2', '<f8', ()),
                                   ('STNMSGM2', '<f8', ()), ('PLINDEX2', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CQUAD8(object):
    name = 'CQUAD8'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4'
    subtables = []


@register_table
class CQUAD8_CN(object):
    name = 'CQUAD8_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4_CN'
    subtables = []


@register_table
class CQUADR(object):
    name = 'CQUADR'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4'
    subtables = []


@register_table
class CQUADR_CN(object):
    name = 'CQUADR_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4_CN'
    subtables = []


@register_table
class CTETRA(object):
    name = 'CTETRA'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CHEXA'
    subtables = []


@register_table
class CTETRA_CN(object):
    name = 'CTETRA_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M0', '<f8', ()), ('M1', '<f8', ()),
                                   ('M2', '<f8', ()), ('M4', '<f8', ()), ('E0', '<f8', ()), ('EP', '<f8', ()),
                                   ('GAMMA', '<f8', ()), ('RMSTRS', '<f8', ()), ('RMSTRN', '<f8', ()),
                                   ('MNSTRS', '<f8', ()), ('MNSTRN', '<f8', ()), ('STRPSGM', '<f8', ()),
                                   ('STRMSGM', '<f8', ()), ('STNPSGM', '<f8', ()), ('STNMSGM', '<f8', ()),
                                   ('PLINDEX', '<f8', ())], (4,)), ('IDENT_ID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIA3(object):
    name = 'CTRIA3'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4'
    subtables = []


@register_table
class CTRIA6(object):
    name = 'CTRIA6'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4'
    subtables = []


@register_table
class CTRIA6_CN(object):
    name = 'CTRIA6_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('EID', '<i8', ()), ('SS',
                                  [('GRID', '<i8', ()), ('LAYER1', '<i8', ()), ('M01', '<f8', ()), ('M11', '<f8', ()),
                                   ('M21', '<f8', ()), ('M41', '<f8', ()), ('E01', '<f8', ()), ('EP1', '<f8', ()),
                                   ('GAMMA1', '<f8', ()), ('RMSTRS1', '<f8', ()), ('RMSTRN1', '<f8', ()),
                                   ('MNSTRS1', '<f8', ()), ('MNSTRN1', '<f8', ()), ('STRPSGM1', '<f8', ()),
                                   ('STRMSGM1', '<f8', ()), ('STNPSGM1', '<f8', ()), ('STNMSGM1', '<f8', ()),
                                   ('PLINDEX1', '<f8', ()), ('LAYER2', '<i8', ()), ('M02', '<f8', ()),
                                   ('M12', '<f8', ()), ('M22', '<f8', ()), ('M42', '<f8', ()), ('E02', '<f8', ()),
                                   ('EP2', '<f8', ()), ('GAMMA2', '<f8', ()), ('RMSTRS2', '<f8', ()),
                                   ('RMSTRN2', '<f8', ()), ('MNSTRS2', '<f8', ()), ('MNSTRN2', '<f8', ()),
                                   ('STRPSGM2', '<f8', ()), ('STRMSGM2', '<f8', ()), ('STNPSGM2', '<f8', ()),
                                   ('STNMSGM2', '<f8', ()), ('PLINDEX2', '<f8', ())], (3,)), ('IDENT_ID', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CTRIAR(object):
    name = 'CTRIAR'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CQUAD4'
    subtables = []


@register_table
class CTRIAR_CN(object):
    name = 'CTRIAR_CN'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/FATIGUE_VIBRATION/CTRIA6_CN'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/FATIGUE_VIBRATION'
    dtype = [('IDENT_ID', '<i8', ()), ('EVNTNAM', 'S56', ()), ('EQNAME', 'S20', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COLUMN(object):
    name = 'COLUMN'
    path = '/NASTRAN/RESULT/MATRIX/BHH'
    dtype = [('POSITION', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/MATRIX/BHH'
    dtype = [('ROW', '<i8', ()), ('VALUE', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/MATRIX/BHH'
    dtype = [('NAME', 'S8', ()), ('FORM', '<i8', ()), ('ROW', '<i8', ()), ('COLUMN', '<i8', ()),
             ('NON_ZERO', '<i8', ()), ('COLUMN_POS', '<i8', ()), ('DATA_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/MATRIX/BHH/COLUMN', '/NASTRAN/RESULT/MATRIX/BHH/DATA']


@register_table
class COLUMN(object):
    name = 'COLUMN'
    path = '/NASTRAN/RESULT/MATRIX/BHH_CPLX'
    dtype = [('POSITION', '<i8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/MATRIX/BHH_CPLX'
    dtype = [('ROW', '<i8', ()), ('RE', '<f8', ()), ('IM', '<f8', ())]
    is_subtable = True
    same_as = 'None'
    subtables = []


@register_table
class IDENTITY(object):
    name = 'IDENTITY'
    path = '/NASTRAN/RESULT/MATRIX/BHH_CPLX'
    dtype = [('NAME', 'S8', ()), ('FORM', '<i8', ()), ('ROW', '<i8', ()), ('COLUMN', '<i8', ()),
             ('NON_ZERO', '<i8', ()), ('COLUMN_POS', '<i8', ()), ('DATA_POS', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = ['/NASTRAN/RESULT/MATRIX/BHH_CPLX/COLUMN', '/NASTRAN/RESULT/MATRIX/BHH_CPLX/DATA']


@register_table
class MCFRAC(object):
    name = 'MCFRAC'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = [('STRMODNO', '<i8', ()), ('STRF_HZ', '<f8', ()), ('MODERESR', '<f8', ()), ('MODERESI', '<f8', ()),
             ('MODERESM', '<f8', ()), ('MODERESP', '<f8', ()), ('PROJMAG', '<f8', ()), ('RELPHASE', '<f8', ()),
             ('MODFRAC', '<f8', ()), ('SRMAG', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SACCELERATION(object):
    name = 'SACCELERATION'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/MODAL/SDISPLACEMENT'
    subtables = []


@register_table
class SACCELERATION_CPLX(object):
    name = 'SACCELERATION_CPLX'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/MODAL/SDISPLACEMENT_CPLX'
    subtables = []


@register_table
class SDISPLACEMENT(object):
    name = 'SDISPLACEMENT'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = [('ID', '<i8', ()), ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SDISPLACEMENT_CPLX(object):
    name = 'SDISPLACEMENT_CPLX'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = [('ID', '<i8', ()), ('RE', '<f8', ()), ('IM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SVECTOR(object):
    name = 'SVECTOR'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NODAL/EIGENVECTOR'
    subtables = []


@register_table
class SVECTOR_CPLX(object):
    name = 'SVECTOR_CPLX'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NODAL/EIGENVECTOR_CPLX'
    subtables = []


@register_table
class SVELOCITY(object):
    name = 'SVELOCITY'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/MODAL/SDISPLACEMENT'
    subtables = []


@register_table
class SVELOCITY_CPLX(object):
    name = 'SVELOCITY_CPLX'
    path = '/NASTRAN/RESULT/MODAL'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/MODAL/SDISPLACEMENT_CPLX'
    subtables = []


@register_table
class AXISYM_1(object):
    name = 'AXISYM_1'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class AXISYM_87(object):
    name = 'AXISYM_87'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class AXISYM_88(object):
    name = 'AXISYM_88'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class AXISYM_89(object):
    name = 'AXISYM_89'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (2,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_14(object):
    name = 'BEAM_14'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class BEAM_52(object):
    name = 'BEAM_52'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class BEAM_78(object):
    name = 'BEAM_78'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class BEAM_79(object):
    name = 'BEAM_79'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class BEAM_98(object):
    name = 'BEAM_98'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class HEXA20_21(object):
    name = 'HEXA20_21'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (27,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA20_35(object):
    name = 'HEXA20_35'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/HEXA20_21'
    subtables = []


@register_table
class HEXA20_57(object):
    name = 'HEXA20_57'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/HEXA8_7'
    subtables = []


@register_table
class HEXA20_61(object):
    name = 'HEXA20_61'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/HEXA8_7'
    subtables = []


@register_table
class HEXA8_117(object):
    name = 'HEXA8_117'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class HEXA8_149(object):
    name = 'HEXA8_149'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class HEXA8_7(object):
    name = 'HEXA8_7'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (8,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFHEXA_188(object):
    name = 'IFHEXA_188'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class IFHEXA_189(object):
    name = 'IFHEXA_189'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class IFPENT_192(object):
    name = 'IFPENT_192'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class IFPENT_193(object):
    name = 'IFPENT_193'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (7,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFQDX_190(object):
    name = 'IFQDX_190'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class IFQDX_191(object):
    name = 'IFQDX_191'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class IFQUAD_186(object):
    name = 'IFQUAD_186'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class IFQUAD_187(object):
    name = 'IFQUAD_187'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class PENTA15_202(object):
    name = 'PENTA15_202'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (21,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA6_136(object):
    name = 'PENTA6_136'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (6,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_10(object):
    name = 'QUAD4_10'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_11(object):
    name = 'QUAD4_11'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_114(object):
    name = 'QUAD4_114'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_115(object):
    name = 'QUAD4_115'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_116(object):
    name = 'QUAD4_116'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_139(object):
    name = 'QUAD4_139'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_140(object):
    name = 'QUAD4_140'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_151(object):
    name = 'QUAD4_151'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class QUAD4_18(object):
    name = 'QUAD4_18'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (4,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_20(object):
    name = 'QUAD4_20'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_3(object):
    name = 'QUAD4_3'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_75(object):
    name = 'QUAD4_75'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_153(object):
    name = 'QUAD8_153'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/AXISYM_89'
    subtables = []


@register_table
class QUAD8_22(object):
    name = 'QUAD8_22'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_26(object):
    name = 'QUAD8_26'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_27(object):
    name = 'QUAD8_27'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_28(object):
    name = 'QUAD8_28'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_30(object):
    name = 'QUAD8_30'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (9,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8_32(object):
    name = 'QUAD8_32'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_33(object):
    name = 'QUAD8_33'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_53(object):
    name = 'QUAD8_53'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_54(object):
    name = 'QUAD8_54'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_55(object):
    name = 'QUAD8_55'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_58(object):
    name = 'QUAD8_58'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_59(object):
    name = 'QUAD8_59'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class QUADX_67(object):
    name = 'QUADX_67'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD8_30'
    subtables = []


@register_table
class ROD_9(object):
    name = 'ROD_9'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class TETRA10_127(object):
    name = 'TETRA10_127'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class TETRA10_130(object):
    name = 'TETRA10_130'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class TETRA10_184(object):
    name = 'TETRA10_184'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class TETRA4_134(object):
    name = 'TETRA4_134'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class TETRA4_157(object):
    name = 'TETRA4_157'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_18'
    subtables = []


@register_table
class TRIA3_138(object):
    name = 'TRIA3_138'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('VALUE', '<f8', (3,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_155(object):
    name = 'TRIA3_155'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA3_156(object):
    name = 'TRIA3_156'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA3_158(object):
    name = 'TRIA3_158'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class TRIA3_2(object):
    name = 'TRIA3_2'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class TRIA3_6(object):
    name = 'TRIA3_6'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/QUAD4_140'
    subtables = []


@register_table
class TRIA6_124(object):
    name = 'TRIA6_124'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_125(object):
    name = 'TRIA6_125'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_126(object):
    name = 'TRIA6_126'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_128(object):
    name = 'TRIA6_128'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_129(object):
    name = 'TRIA6_129'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_200(object):
    name = 'TRIA6_200'
    path = '/NASTRAN/RESULT/NLOUT/SCALAR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/SCALAR/IFPENT_193'
    subtables = []


@register_table
class AXISYM_1(object):
    name = 'AXISYM_1'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class AXISYM_87(object):
    name = 'AXISYM_87'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class AXISYM_88(object):
    name = 'AXISYM_88'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class AXISYM_89(object):
    name = 'AXISYM_89'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (12,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class BEAM_14(object):
    name = 'BEAM_14'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class BEAM_52(object):
    name = 'BEAM_52'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class BEAM_78(object):
    name = 'BEAM_78'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class BEAM_79(object):
    name = 'BEAM_79'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class BEAM_98(object):
    name = 'BEAM_98'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class HEXA20_21(object):
    name = 'HEXA20_21'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (162,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HEXA20_35(object):
    name = 'HEXA20_35'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/HEXA20_21'
    subtables = []


@register_table
class HEXA20_57(object):
    name = 'HEXA20_57'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/HEXA8_7'
    subtables = []


@register_table
class HEXA20_61(object):
    name = 'HEXA20_61'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/HEXA8_7'
    subtables = []


@register_table
class HEXA8_117(object):
    name = 'HEXA8_117'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class HEXA8_149(object):
    name = 'HEXA8_149'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class HEXA8_7(object):
    name = 'HEXA8_7'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (48,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFHEXA_188(object):
    name = 'IFHEXA_188'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class IFHEXA_189(object):
    name = 'IFHEXA_189'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class IFPENT_192(object):
    name = 'IFPENT_192'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class IFPENT_193(object):
    name = 'IFPENT_193'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (42,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IFQDX_190(object):
    name = 'IFQDX_190'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class IFQDX_191(object):
    name = 'IFQDX_191'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class IFQUAD_186(object):
    name = 'IFQUAD_186'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class IFQUAD_187(object):
    name = 'IFQUAD_187'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class PENTA15_202(object):
    name = 'PENTA15_202'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (126,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PENTA6_136(object):
    name = 'PENTA6_136'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (36,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_10(object):
    name = 'QUAD4_10'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_11(object):
    name = 'QUAD4_11'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_114(object):
    name = 'QUAD4_114'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_115(object):
    name = 'QUAD4_115'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_116(object):
    name = 'QUAD4_116'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class QUAD4_139(object):
    name = 'QUAD4_139'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_140(object):
    name = 'QUAD4_140'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (6,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_151(object):
    name = 'QUAD4_151'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class QUAD4_18(object):
    name = 'QUAD4_18'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (24,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD4_20(object):
    name = 'QUAD4_20'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_3(object):
    name = 'QUAD4_3'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD4_75(object):
    name = 'QUAD4_75'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_153(object):
    name = 'QUAD8_153'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/AXISYM_89'
    subtables = []


@register_table
class QUAD8_22(object):
    name = 'QUAD8_22'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_26(object):
    name = 'QUAD8_26'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_27(object):
    name = 'QUAD8_27'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_28(object):
    name = 'QUAD8_28'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_30(object):
    name = 'QUAD8_30'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (54,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class QUAD8_32(object):
    name = 'QUAD8_32'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_33(object):
    name = 'QUAD8_33'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class QUAD8_53(object):
    name = 'QUAD8_53'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_54(object):
    name = 'QUAD8_54'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_55(object):
    name = 'QUAD8_55'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_58(object):
    name = 'QUAD8_58'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUAD8_59(object):
    name = 'QUAD8_59'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class QUADX_67(object):
    name = 'QUADX_67'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD8_30'
    subtables = []


@register_table
class ROD_9(object):
    name = 'ROD_9'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class TETRA10_127(object):
    name = 'TETRA10_127'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class TETRA10_130(object):
    name = 'TETRA10_130'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class TETRA10_184(object):
    name = 'TETRA10_184'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class TETRA4_134(object):
    name = 'TETRA4_134'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class TETRA4_157(object):
    name = 'TETRA4_157'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_18'
    subtables = []


@register_table
class TRIA3_138(object):
    name = 'TRIA3_138'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = [('EID', '<i8', ()), ('PLYID', '<i8', ()), ('POSTID', '<i8', ()), ('POSTNAME', 'S8', ()),
             ('CID', '<i8', ()), ('VALUE', '<f8', (18,)), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TRIA3_155(object):
    name = 'TRIA3_155'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA3_156(object):
    name = 'TRIA3_156'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA3_158(object):
    name = 'TRIA3_158'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class TRIA3_2(object):
    name = 'TRIA3_2'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class TRIA3_6(object):
    name = 'TRIA3_6'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/QUAD4_140'
    subtables = []


@register_table
class TRIA6_124(object):
    name = 'TRIA6_124'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_125(object):
    name = 'TRIA6_125'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_126(object):
    name = 'TRIA6_126'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_128(object):
    name = 'TRIA6_128'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_129(object):
    name = 'TRIA6_129'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/TRIA3_138'
    subtables = []


@register_table
class TRIA6_200(object):
    name = 'TRIA6_200'
    path = '/NASTRAN/RESULT/NLOUT/TENSOR'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NLOUT/TENSOR/IFPENT_193'
    subtables = []


@register_table
class ACCELERATION(object):
    name = 'ACCELERATION'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACCELERATION_CPLX(object):
    name = 'ACCELERATION_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACCELERATION_TRANS(object):
    name = 'ACCELERATION_TRANS'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACCELERATION_TRANS_CPLX(object):
    name = 'ACCELERATION_TRANS_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('XI', '<f8', ()),
             ('YI', '<f8', ()), ('ZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class APPLIED_LOAD(object):
    name = 'APPLIED_LOAD'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class APPLIED_LOAD_CPLX(object):
    name = 'APPLIED_LOAD_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DISPLACEMENT(object):
    name = 'DISPLACEMENT'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DISPLACEMENT_CPLX(object):
    name = 'DISPLACEMENT_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DISPLACEMENT_TRANS(object):
    name = 'DISPLACEMENT_TRANS'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DISPLACEMENT_TRANS_CPLX(object):
    name = 'DISPLACEMENT_TRANS_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('XI', '<f8', ()),
             ('YI', '<f8', ()), ('ZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVECTOR(object):
    name = 'EIGENVECTOR'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVECTOR_CPLX(object):
    name = 'EIGENVECTOR_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVECTOR_TRANS(object):
    name = 'EIGENVECTOR_TRANS'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVECTOR_TRANS_CPLX(object):
    name = 'EIGENVECTOR_TRANS_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('XI', '<f8', ()),
             ('YI', '<f8', ()), ('ZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_FORCE(object):
    name = 'GRID_FORCE'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('EID', '<i8', ()), ('ELNAME', 'S8', ()), ('F1', '<f8', ()), ('F2', '<f8', ()),
             ('F3', '<f8', ()), ('M1', '<f8', ()), ('M2', '<f8', ()), ('M3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_FORCE_CPLX(object):
    name = 'GRID_FORCE_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('EID', '<i8', ()), ('ELNAME', 'S8', ()), ('F1R', '<f8', ()), ('F2R', '<f8', ()),
             ('F3R', '<f8', ()), ('M1R', '<f8', ()), ('M2R', '<f8', ()), ('M3R', '<f8', ()), ('F1I', '<f8', ()),
             ('F2I', '<f8', ()), ('F3I', '<f8', ()), ('M1I', '<f8', ()), ('M2I', '<f8', ()), ('M3I', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_WEIGHT(object):
    name = 'GRID_WEIGHT'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('MO', '<f8', (36,)), ('S', '<f8', (9,)), ('MX', '<f8', ()), ('XX', '<f8', ()),
             ('YX', '<f8', ()), ('ZX', '<f8', ()), ('MY', '<f8', ()), ('XY', '<f8', ()), ('YY', '<f8', ()),
             ('ZY', '<f8', ()), ('MZ', '<f8', ()), ('XZ', '<f8', ()), ('YZ', '<f8', ()), ('ZZ', '<f8', ()),
             ('I', '<f8', (9,)), ('PIX', '<f8', ()), ('PIY', '<f8', ()), ('PIZ', '<f8', ()), ('Q', '<f8', (9,)),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class KINETIC_ENERGY(object):
    name = 'KINETIC_ENERGY'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('TYPE', '<i8', ()), ('KET1', '<f8', ()), ('KET2', '<f8', ()), ('KET3', '<f8', ()),
             ('KER1', '<f8', ()), ('KER2', '<f8', ()), ('KER3', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class KINETIC_ENERGY_CPLX(object):
    name = 'KINETIC_ENERGY_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('TYPE', '<i8', ()), ('KET1R', '<f8', ()), ('KET2R', '<f8', ()), ('KET3R', '<f8', ()),
             ('KER1R', '<f8', ()), ('KER2R', '<f8', ()), ('KER3R', '<f8', ()), ('KET1I', '<f8', ()),
             ('KET2I', '<f8', ()), ('KET3I', '<f8', ()), ('KER1I', '<f8', ()), ('KER2I', '<f8', ()),
             ('KER3I', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MPC_FORCE(object):
    name = 'MPC_FORCE'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MPC_FORCE_CPLX(object):
    name = 'MPC_FORCE_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACCELERATION_CPLX(object):
    name = 'ACCELERATION_CPLX'
    path = '/NASTRAN/RESULT/NODAL/RESPONSE_SPECTRUM'
    dtype = [('ID', '<i8', ()), ('FREQUNCY', '<f8', ()), ('T1R', '<f8', ()), ('T2R', '<f8', ()), ('T3R', '<f8', ()),
             ('R1R', '<f8', ()), ('R2R', '<f8', ()), ('R3R', '<f8', ()), ('T1I', '<f8', ()), ('T2I', '<f8', ()),
             ('T3I', '<f8', ()), ('R1I', '<f8', ()), ('R2I', '<f8', ()), ('R3I', '<f8', ()), ('DAMP', '<f8', ()),
             ('RECORD', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DISPLACEMENT_CPLX(object):
    name = 'DISPLACEMENT_CPLX'
    path = '/NASTRAN/RESULT/NODAL/RESPONSE_SPECTRUM'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NODAL/RESPONSE_SPECTRUM/ACCELERATION_CPLX'
    subtables = []


@register_table
class VELOCITY_CPLX(object):
    name = 'VELOCITY_CPLX'
    path = '/NASTRAN/RESULT/NODAL/RESPONSE_SPECTRUM'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NODAL/RESPONSE_SPECTRUM/ACCELERATION_CPLX'
    subtables = []


@register_table
class SQ_WETTED(object):
    name = 'SQ_WETTED'
    path = '/NASTRAN/RESULT/NODAL/SENSITIVITY'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/NODAL/SENSITIVITY/WETTED'
    subtables = []


@register_table
class WETTED(object):
    name = 'WETTED'
    path = '/NASTRAN/RESULT/NODAL/SENSITIVITY'
    dtype = [('ID', '<i8', ()), ('RESPDOF', '<i8', ()), ('COMP', '<i8', ()), ('SENSR', '<f8', ()), ('SENSI', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPC_FORCE(object):
    name = 'SPC_FORCE'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SPC_FORCE_CPLX(object):
    name = 'SPC_FORCE_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELEM_DISCD_SURF(object):
    name = 'ELEM_DISCD_SURF'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('ELTYPE', 'S8', ()), ('FIBRE', 'S4', ()), ('NX', '<f8', ()), ('NY', '<f8', ()),
             ('TXY', '<f8', ()), ('MJRP', '<f8', ()), ('MNRP', '<f8', ()), ('TMAX', '<f8', ()), ('HVM', '<f8', ()),
             ('ERR', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELEM_DISCD_VOL_DCT(object):
    name = 'ELEM_DISCD_VOL_DCT'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('ELTYPE', 'S8', ()), ('NX', '<f8', ()), ('NY', '<f8', ()), ('NZ', '<f8', ()),
             ('SXY', '<f8', ()), ('SYZ', '<f8', ()), ('SZX', '<f8', ()), ('MP', '<f8', ()), ('HVM', '<f8', ()),
             ('ERR', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ELEM_DISCD_VOL_PRIN(object):
    name = 'ELEM_DISCD_VOL_PRIN'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('ELTYPE', 'S8', ()), ('SA', '<f8', ()), ('SB', '<f8', ()), ('SC', '<f8', ()),
             ('MP', '<f8', ()), ('HVM', '<f8', ()), ('ERR', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_DISCD_SURF(object):
    name = 'GRID_DISCD_SURF'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('FIBRE', 'S4', ()), ('NX', '<f8', ()), ('NY', '<f8', ()), ('NXY', '<f8', ()),
             ('MJRP', '<f8', ()), ('MNRP', '<f8', ()), ('TMAX', '<f8', ()), ('HVM', '<f8', ()), ('ERR', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_DISCD_VOL_DCT(object):
    name = 'GRID_DISCD_VOL_DCT'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('NX', '<f8', ()), ('NY', '<f8', ()), ('NZ', '<f8', ()), ('TXY', '<f8', ()),
             ('TYZ', '<f8', ()), ('TZX', '<f8', ()), ('PR', '<f8', ()), ('HVM', '<f8', ()), ('ERR', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_DISCD_VOL_PRIN(object):
    name = 'GRID_DISCD_VOL_PRIN'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('SA', '<f8', ()), ('SB', '<f8', ()), ('SC', '<f8', ()), ('PR', '<f8', ()),
             ('HVM', '<f8', ()), ('ERR', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRID_SURF_PLANE_STRAIN(object):
    name = 'GRID_SURF_PLANE_STRAIN'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('NX', '<f8', ()), ('NY', '<f8', ()), ('NZ', '<f8', ()), ('TXY', '<f8', ()),
             ('PR', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('IDENT', '<i8', ()), ('ID', '<i8', ()), ('REFID', '<i8', ()), ('OCOORD', '<i8', ()), ('AXIS', '<i8', ()),
             ('NORMAL', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class SURFACE(object):
    name = 'SURFACE'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('EID', '<i8', ()), ('FIBRE', 'S4', ()), ('NX', '<f8', ()), ('NY', '<f8', ()),
             ('TXY', '<f8', ()), ('A', '<f8', ()), ('MJRP', '<f8', ()), ('MNRP', '<f8', ()), ('TMAX', '<f8', ()),
             ('HVM', '<f8', ()), ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VOLUME_DIRECT(object):
    name = 'VOLUME_DIRECT'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('NX', '<f8', ()), ('NY', '<f8', ()), ('NZ', '<f8', ()), ('TXY', '<f8', ()),
             ('TYZ', '<f8', ()), ('TZX', '<f8', ()), ('PR', '<f8', ()), ('HVM', '<f8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VOLUME_PRINCIPAL(object):
    name = 'VOLUME_PRINCIPAL'
    path = '/NASTRAN/RESULT/NODAL/STRESS'
    dtype = [('ID', '<i8', ()), ('LXA', '<f8', ()), ('LXB', '<f8', ()), ('LXC', '<f8', ()), ('LYA', '<f8', ()),
             ('LYB', '<f8', ()), ('LYC', '<f8', ()), ('LZA', '<f8', ()), ('LZB', '<f8', ()), ('LZC', '<f8', ()),
             ('SA', '<f8', ()), ('SB', '<f8', ()), ('SC', '<f8', ()), ('EPR', '<f8', ()), ('EHVM', '<f8', ()),
             ('IDENT', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPERATURE(object):
    name = 'TEMPERATURE'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TEMPERATURE_SHELL(object):
    name = 'TEMPERATURE_SHELL'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('T_TOP', '<f8', ()), ('T_BOT', '<f8', ()), ('T_MID', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VELOCITY(object):
    name = 'VELOCITY'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('RX', '<f8', ()),
             ('RY', '<f8', ()), ('RZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VELOCITY_CPLX(object):
    name = 'VELOCITY_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('RXR', '<f8', ()),
             ('RYR', '<f8', ()), ('RZR', '<f8', ()), ('XI', '<f8', ()), ('YI', '<f8', ()), ('ZI', '<f8', ()),
             ('RXI', '<f8', ()), ('RYI', '<f8', ()), ('RZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VELOCITY_TRANS(object):
    name = 'VELOCITY_TRANS'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('X', '<f8', ()), ('Y', '<f8', ()), ('Z', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VELOCITY_TRANS_CPLX(object):
    name = 'VELOCITY_TRANS_CPLX'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('ID', '<i8', ()), ('XR', '<f8', ()), ('YR', '<f8', ()), ('ZR', '<f8', ()), ('XI', '<f8', ()),
             ('YI', '<f8', ()), ('ZI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VELOCITY_NORMAL(object):
    name = 'VELOCITY_NORMAL'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('GID', '<i8', ()), ('UNX', '<f8', ()), ('UNY', '<f8', ()), ('UNZ', '<f8', ()), ('NVALR', '<f8', ()),
             ('NVALI', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class VIRTUAL_CRACK(object):
    name = 'VIRTUAL_CRACK'
    path = '/NASTRAN/RESULT/NODAL'
    dtype = [('CID', '<i8', ()), ('GID', '<i8', ()), ('ERATE_TOTAL', '<f8', ()), ('ERATE1', '<f8', ()),
             ('ERATE2', '<f8', ()), ('ERATE3', '<f8', ()), ('DX', '<f8', ()), ('DY', '<f8', ()), ('DZ', '<f8', ()),
             ('XX', '<f8', ()), ('XY', '<f8', ()), ('XZ', '<f8', ()), ('YX', '<f8', ()), ('YY', '<f8', ()),
             ('YZ', '<f8', ()), ('ZX', '<f8', ()), ('ZY', '<f8', ()), ('ZZ', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CONSTRAINT(object):
    name = 'CONSTRAINT'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('IDCID', '<i8', ()), ('DCID', '<i8', ()), ('IRID', '<i8', ()), ('RTYPE', '<i8', ()), ('TYPE', '<i8', ()),
             ('LUFLAG', '<i8', ()), ('BOUND', '<f8', ()), ('REGION', '<i8', ()), ('VALUE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GRIDNEW(object):
    name = 'GRIDNEW'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('ID', '<i8', ()), ('CP', '<i8', ()), ('X', '<f8', (3,)), ('CD', '<i8', ()), ('PS', '<i8', ()),
             ('SEID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class HISTORY(object):
    name = 'HISTORY'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('VID', '<i8', ()), ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LABEL(object):
    name = 'LABEL'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('IDVID', '<i8', ()), ('DVID', '<i8', ()), ('LABEL', 'S8', ()), ('VMIN', '<f8', ()), ('VMAX', '<f8', ()),
             ('DELX', '<f8', ()), ('DVSID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class OBJECTIVE(object):
    name = 'OBJECTIVE'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('EXACT', '<f8', ()), ('APPRX', '<f8', ()), ('MAXIM', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ABSTRESS(object):
    name = 'ABSTRESS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class ACINTS(object):
    name = 'ACINTS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('WGID', '<i8', ()), ('RESPONSE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ACPWR(object):
    name = 'ACPWR'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('PANELNAME', 'S8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AFPRES(object):
    name = 'AFPRES'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('AFPCOMP', '<i8', ()),
             ('AFPMGID', '<i8', ()), ('AFPMID', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AFPWR(object):
    name = 'AFPWR'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('AFPMID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class AFINTS(object):
    name = 'AFINTS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/AFPRES'
    subtables = []


@register_table
class AFVELO(object):
    name = 'AFVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/AFPRES'
    subtables = []


@register_table
class CEIG(object):
    name = 'CEIG'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('AORO', '<i8', ()), ('RESPONSE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CFAILURE(object):
    name = 'CFAILURE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ICODE', '<i8', ()), ('PLY', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COMP(object):
    name = 'COMP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/VOLUME'
    subtables = []


@register_table
class CSTRAT(object):
    name = 'CSTRAT'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/CFAILURE'
    subtables = []


@register_table
class CSTRAIN(object):
    name = 'CSTRAIN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/CFAILURE'
    subtables = []


@register_table
class CSTRESS(object):
    name = 'CSTRESS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/CFAILURE'
    subtables = []


@register_table
class DISP(object):
    name = 'DISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('COMP', '<i8', ()), ('GRID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DIVERG(object):
    name = 'DIVERG'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ROOT', '<i8', ()), ('MACH', '<f8', ()),
             ('DIVERG', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGN(object):
    name = 'EIGN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('APRX', '<i8', ()), ('STRUFLUD', 'S8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ERP(object):
    name = 'ERP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ICODE', '<i8', ()), ('SET3ID', '<i8', ()),
             ('PANELNAME', 'S8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class ESE(object):
    name = 'ESE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ESECOMP', '<i8', ()), ('ESEID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FATIGUE(object):
    name = 'FATIGUE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ICODE', '<i8', ()), ('FTGID', '<i8', ()),
             ('EID', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLUTTER(object):
    name = 'FLUTTER'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('MACH', '<f8', ()), ('VELOCITY', '<f8', ()),
             ('DENSITY', '<f8', ()), ('FLUTTER', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCE(object):
    name = 'FORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class FRFTG(object):
    name = 'FRFTG'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/FATIGUE'
    subtables = []


@register_table
class FRACCL(object):
    name = 'FRACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class FRDISP(object):
    name = 'FRDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class FRFORC(object):
    name = 'FRFORC'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class FRSPCF(object):
    name = 'FRSPCF'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class FRSTRE(object):
    name = 'FRSTRE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class FRVELO(object):
    name = 'FRVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class FREQ(object):
    name = 'FREQ'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/EIGN'
    subtables = []


@register_table
class FRMASS(object):
    name = 'FRMASS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('TOPVARID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GMASS(object):
    name = 'GMASS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/EIGN'
    subtables = []


@register_table
class GPFORCE(object):
    name = 'GPFORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('GPFCOMP', '<i8', ()),
             ('GPFEID', '<i8', ()), ('GPFGID', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GPFORCP(object):
    name = 'GPFORCP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('GPRGID', '<i8', ()), ('GPPGID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class GSTIFF(object):
    name = 'GSTIFF'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/EIGN'
    subtables = []


@register_table
class LAMA(object):
    name = 'LAMA'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('APRX', '<i8', ()), ('RESPONSE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PRES(object):
    name = 'PRES'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ACCOMP', '<i8', ()), ('ACGRID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSDDISP(object):
    name = 'PSDDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('COMP', '<i8', ()), ('RANDPS', '<i8', ()),
             ('GRID', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSDVELO(object):
    name = 'PSDVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/PSDDISP'
    subtables = []


@register_table
class PSDACCL(object):
    name = 'PSDACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSDISP(object):
    name = 'RMSDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSVELO(object):
    name = 'RMSVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSACCL(object):
    name = 'RMSACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/PSDDISP'
    subtables = []


@register_table
class SPCFORCE(object):
    name = 'SPCFORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class STABDER(object):
    name = 'STABDER'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('XID', '<i8', ()), ('COMP', '<i8', ()),
             ('RESFLG', '<i8', ()), ('TRIM', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class STRAIN(object):
    name = 'STRAIN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class STRESS(object):
    name = 'STRESS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('ICODE', '<i8', ()), ('EPID', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TACCL(object):
    name = 'TACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class TDISP(object):
    name = 'TDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class TFORC(object):
    name = 'TFORC'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class TOTSE(object):
    name = 'TOTSE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/VOLUME'
    subtables = []


@register_table
class TRIM(object):
    name = 'TRIM'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('XID', '<i8', ()), ('TRIM', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TSPCF(object):
    name = 'TSPCF'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class TSTRE(object):
    name = 'TSTRE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/STRESS'
    subtables = []


@register_table
class TVELO(object):
    name = 'TVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1/DISP'
    subtables = []


@register_table
class VOLUME(object):
    name = 'VOLUME'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('RESPONSE', '<f8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WEIGHT(object):
    name = 'WEIGHT'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('RBWROW', '<i8', ()), ('RBWCOL', '<i8', ()),
             ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WMPID(object):
    name = 'WMPID'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE/RTYPE1'
    dtype = [('IRID', '<i8', ()), ('RID', '<i8', ()), ('LABEL', 'S8', ()), ('MATID', '<i8', ()), ('PROPID', '<i8', ()),
             ('PTYPE', '<i8', ()), ('WMPRNG', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTYPE2(object):
    name = 'RTYPE2'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE'
    dtype = [('IR2ID', '<i8', ()), ('R2ID', '<i8', ()), ('LABEL', 'S8', ()), ('EQID', '<i8', ()), ('REGION', '<i8', ()),
             ('METHOD', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTYPE3(object):
    name = 'RTYPE3'
    path = '/NASTRAN/RESULT/OPTIMIZATION/RESPONSE'
    dtype = [('IR3ID', '<i8', ()), ('R3ID', '<i8', ()), ('LABEL', 'S8', ()), ('GROUP', 'S8', ()), ('TYPE', 'S8', ()),
             ('REGION', '<i8', ()), ('RESPONSE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COEFFICIENT(object):
    name = 'COEFFICIENT'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY'
    dtype = [('VID', '<i8', ()), ('COEFFICIENT', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class CEIG(object):
    name = 'CEIG'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class CFAILURE(object):
    name = 'CFAILURE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('EID', '<i8', ()), ('COMP', '<i8', ()), ('PLY', '<i8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class COMP(object):
    name = 'COMP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class CSTRAIN(object):
    name = 'CSTRAIN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/CFAILURE'
    subtables = []


@register_table
class CSTRESS(object):
    name = 'CSTRESS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/CFAILURE'
    subtables = []


@register_table
class DISP(object):
    name = 'DISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('GRID', '<i8', ()), ('COMP', '<i8', ()), ('POS', '<i8', ()),
             ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DIVERG(object):
    name = 'DIVERG'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('ROOT', '<i8', ()), ('MACH', '<f8', ()), ('POS', '<i8', ()),
             ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGN(object):
    name = 'EIGN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class ESE(object):
    name = 'ESE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('EID', '<i8', ()), ('COMP', '<i8', ()), ('POS', '<i8', ()),
             ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FLUTTER(object):
    name = 'FLUTTER'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('DENSITY', '<f8', ()), ('MACH', '<f8', ()),
             ('VELOCITY', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FORCE(object):
    name = 'FORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/STRESS'
    subtables = []


@register_table
class FRACCL(object):
    name = 'FRACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class FRDISP(object):
    name = 'FRDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class FREQ(object):
    name = 'FREQ'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class FRFORC(object):
    name = 'FRFORC'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/ESE'
    subtables = []


@register_table
class FRMASS(object):
    name = 'FRMASS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('TOPVAR', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class FRSPCF(object):
    name = 'FRSPCF'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class FRSTRE(object):
    name = 'FRSTRE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/ESE'
    subtables = []


@register_table
class FRVELO(object):
    name = 'FRVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class GPFORCE(object):
    name = 'GPFORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/ESE'
    subtables = []


@register_table
class GPFORCP(object):
    name = 'GPFORCP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('GRID', '<i8', ()), ('ORIENT', '<i8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class LAMA(object):
    name = 'LAMA'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class PRES(object):
    name = 'PRES'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class PSDACCL(object):
    name = 'PSDACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/PSDDISP'
    subtables = []


@register_table
class PSDDISP(object):
    name = 'PSDDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('GRID', '<i8', ()), ('COMP', '<i8', ()),
             ('RANDPS', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class PSDVELO(object):
    name = 'PSDVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSACCL(object):
    name = 'RMSACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSDISP(object):
    name = 'RMSDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/PSDDISP'
    subtables = []


@register_table
class RMSVELO(object):
    name = 'RMSVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/PSDDISP'
    subtables = []


@register_table
class SPCFORCE(object):
    name = 'SPCFORCE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class STABDER(object):
    name = 'STABDER'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('RU', '<i8', ()), ('COMP', '<i8', ()), ('POS', '<i8', ()),
             ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class STRAIN(object):
    name = 'STRAIN'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/STRESS'
    subtables = []


@register_table
class STRESS(object):
    name = 'STRESS'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('EID', '<i8', ()), ('COMP', '<i8', ()),
             ('VIEWID', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TACCL(object):
    name = 'TACCL'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class TDISP(object):
    name = 'TDISP'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class TFORC(object):
    name = 'TFORC'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/ESE'
    subtables = []


@register_table
class TOTSE(object):
    name = 'TOTSE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class TSPCF(object):
    name = 'TSPCF'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class TSTRE(object):
    name = 'TSTRE'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/ESE'
    subtables = []


@register_table
class TRIM(object):
    name = 'TRIM'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('XID', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TVELO(object):
    name = 'TVELO'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/DISP'
    subtables = []


@register_table
class VOLUME(object):
    name = 'VOLUME'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1/WEIGHT'
    subtables = []


@register_table
class WEIGHT(object):
    name = 'WEIGHT'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class WMPID(object):
    name = 'WMPID'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY/RTYPE1'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('MID', '<i8', ()), ('PID', '<i8', ()), ('PTYPE', '<i8', ()),
             ('POS', '<i8', ()), ('LEN', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTYPE2(object):
    name = 'RTYPE2'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class RTYPE3(object):
    name = 'RTYPE3'
    path = '/NASTRAN/RESULT/OPTIMIZATION/SENSITIVITY'
    dtype = [('RID', '<i8', ()), ('RESPONSE', '<f8', ()), ('RESPN', '<i8', ()), ('POS', '<i8', ()), ('LEN', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TOMVAR(object):
    name = 'TOMVAR'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('EID', '<i8', ()), ('VALUE', '<f8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class TOPVAR(object):
    name = 'TOPVAR'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = []
    is_subtable = False
    same_as = '/NASTRAN/RESULT/OPTIMIZATION/TOMVAR'
    subtables = []


@register_table
class VARIABLE(object):
    name = 'VARIABLE'
    path = '/NASTRAN/RESULT/OPTIMIZATION'
    dtype = [('VID', '<i8', ()), ('INIT', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVALUE(object):
    name = 'EIGENVALUE'
    path = '/NASTRAN/RESULT/SUMMARY'
    dtype = [('MODE', '<i8', ()), ('ORDER', '<i8', ()), ('EIGEN', '<f8', ()), ('OMEGA', '<f8', ()), ('FREQ', '<f8', ()),
             ('MASS', '<f8', ()), ('STIFF', '<f8', ()), ('RESFLG', '<i8', ()), ('FLDFLG', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class EIGENVALUE_CPLX(object):
    name = 'EIGENVALUE_CPLX'
    path = '/NASTRAN/RESULT/SUMMARY'
    dtype = [('MODE', '<i8', ()), ('ORDER', '<i8', ()), ('REIGEN', '<f8', ()), ('IEIGEN', '<f8', ()),
             ('FREQ', '<f8', ()), ('DAMP', '<f8', ()), ('SPIN_SPEED', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class DATA(object):
    name = 'DATA'
    path = '/NASTRAN/RESULT/SUMMARY/MAX_MIN'
    dtype = [('ID', '<i8', ()), ('CID', '<i8', ()), ('OPER', 'S4', ()), ('DEPTH', '<i8', ()), ('CNAME', 'S8', ()),
             ('IVALUE', '<i8', ()), ('RVALUE', '<f8', ()), ('VALUE', '<f8', ()), ('GRID', '<i8', ()),
             ('SDIST', '<f8', ()), ('RMS', '<f8', ()), ('NRMS', '<i8', ()), ('IDENT', '<i8', ()),
             ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/SUMMARY/MAX_MIN'
    dtype = [('IDENT', '<i8', ()), ('ACODE', '<i8', ()), ('ELTYPE', '<i8', ()), ('ELFLAG', '<i8', ()),
             ('SECAPP', '<i8', ()), ('ASSVAR', '<i8', ()), ('RMSFLG', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class IDENT(object):
    name = 'IDENT'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('ACODE', '<i8', ()), ('SUBCASE', '<i8', ()), ('FRTM', '<f8', ()), ('NUMWDE', '<i8', ()),
             ('CLASS', '<i8', ()), ('NAME', 'S8', ()), ('CP', '<i8', ()), ('CD', '<i8', ()), ('X', '<f8', ()),
             ('Y', '<f8', ()), ('Z', '<f8', ()), ('GRIDID', '<i8', ()), ('ELEMID', '<i8', ()), ('COMP', 'S8', ()),
             ('MESH', '<i8', ()), ('MLABEL', 'S56', ()), ('MPID', '<i8', ()), ('DOMAIN_ID', '<i8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON1_AESTAT_AERO(object):
    name = 'MON1_AESTAT_AERO'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('AXIS', '<i8', ()), ('RIGAIR', '<f8', ()), ('ELARES', '<f8', ()), ('RESAPP', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON1_AESTAT_STRU(object):
    name = 'MON1_AESTAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('AXIS', '<i8', ()), ('RIGAIR', '<f8', ()), ('ELARES', '<f8', ()), ('INER', '<f8', ()),
             ('RIGAPP', '<f8', ()), ('RESAPP', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON1_FREQ_STRU(object):
    name = 'MON1_FREQ_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('AXIS', '<i8', ()), ('INERR', '<f8', ()), ('EXTR', '<f8', ()), ('TOTALR', '<f8', ()),
             ('INERI', '<f8', ()), ('EXTI', '<f8', ()), ('TOTALI', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON1_STAT_STRU(object):
    name = 'MON1_STAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('RIGAPP', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON1_TRAN_STRU(object):
    name = 'MON1_TRAN_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('AXIS', '<i8', ()), ('INER', '<f8', ()), ('EXT', '<f8', ()), ('TOTAL', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON3_AESTAT_STRU(object):
    name = 'MON3_AESTAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('ELARES', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON3_FREQ_STRU(object):
    name = 'MON3_FREQ_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('RESULTR', '<f8', ()), ('RESULTI', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON3_STAT_STRU(object):
    name = 'MON3_STAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('RESAPP', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MON3_TRAN_STRU(object):
    name = 'MON3_TRAN_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('RESULT', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MONA_TRAN_STRU(object):
    name = 'MONA_TRAN_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('A', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOND_AESTAT_AERO(object):
    name = 'MOND_AESTAT_AERO'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('ELARES', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOND_AESTAT_STRU(object):
    name = 'MOND_AESTAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('ELARES', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOND_FREQ_STRU(object):
    name = 'MOND_FREQ_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('DR', '<f8', ()), ('DI', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOND_STAT_STRU(object):
    name = 'MOND_STAT_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('ELARES', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MOND_TRAN_STRU(object):
    name = 'MOND_TRAN_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('D', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []


@register_table
class MONV_TRAN_STRU(object):
    name = 'MONV_TRAN_STRU'
    path = '/NASTRAN/RESULT/MONITOR'
    dtype = [('V', '<f8', ())]
    is_subtable = False
    same_as = 'None'
    subtables = []
