#pylint: disable=R0913
"""
Defines the sub-OP2 class.  This should never be called outside of the OP2 class.

 - OP2_Scalar(debug=False, log=None, debug_file=None)

   **Methods**
   - set_subcases(subcases=None)
   - set_transient_times(times)
   - read_op2(op2_filename=None, combine=False)
   - set_additional_generalized_tables_to_read(tables)
   - set_additional_result_tables_to_read(tables)
   - set_additional_matrices_to_read(matrices)

   **Attributes**
   - total_effective_mass_matrix
   - effective_mass_matrix
   - rigid_body_mass_matrix
   - modal_effective_mass_fraction
   - modal_participation_factors
   - modal_effective_mass
   - modal_effective_weight
   - set_as_msc()
   - set_as_optistruct()

   **Private Methods**
   - _get_table_mapper()
   - _not_available(data, ndata)
   - _table_crasher(data, ndata)
   - _table_passer(data, ndata)
   - _validate_op2_filename(op2_filename)
   - _create_binary_debug()
   - _make_tables()
   - _read_tables(table_name)
   - _skip_table(table_name)
   - _read_table_name(rewind=False, stop_on_failure=True)
   - _update_generalized_tables(tables)
   - _read_cmodext()
   - _read_cmodext_helper(marker_orig, debug=False)
   - _read_geom_table()
   - _finish()

"""
import os
from struct import Struct, unpack
from collections import defaultdict
from typing import List, Tuple, Dict, Union, Optional, Any

from numpy import array
import numpy as np
import cpylog
if cpylog.__version__ >= '1.5.0':  # pragma: no cover
    #import warnings
    #warnings.warn('run "pip install cpylog>=1.5.0"')
    from cpylog import get_logger2, log_exc
else:  # pramga: no cover
    from cpylog import get_logger2
    def log_exc(*args, **kwargs):
        pass


from pyNastran import is_release, __version__
from pyNastran.f06.errors import FatalError
from pyNastran.op2.errors import EmptyRecordError
from pyNastran.op2.op2_interface.op2_reader import OP2Reader, reshape_bytes_block
from pyNastran.bdf.cards.params import PARAM

#============================

from pyNastran.op2.op2_interface.msc_tables import MSC_RESULT_TABLES, MSC_MATRIX_TABLES, MSC_GEOM_TABLES
from pyNastran.op2.op2_interface.nx_tables import NX_RESULT_TABLES, NX_MATRIX_TABLES, NX_GEOM_TABLES

from pyNastran.op2.op2_interface.op2_common import OP2Common
from pyNastran.op2.fortran_format import FortranFormat

from pyNastran.utils import is_binary_file
"""
ftp://161.24.15.247/Nastran2011/seminar/SEC04-DMAP_MODULES.pdf

Datablock	Type	Description
EFMFSMS	Matrix	6 x 1 Total Effective mass matrix
EFMASSS	Matrix	6 x 6 Effective mass matrix
RBMASS	Matrix	6 x 6 Rigid body mass matrix
EFMFACS	Matrix	6 X N Modal effective mass fraction matrix
MPFACS	Matrix	6 x N Modal participation factor matrix
MEFMASS	Matrix	6 x N Modal effective mass matrix
MEFWTS	Matrix	6 x N Modal effective weight matrix
RAFGEN	Matrix	N x M Generalized force matrix
RADEFMP	Matrix	N X U2 Effective inertia loads
BHH	Matrix	N x N Viscous damping matrix
K4HH	Matrix	N x N Structural damping matrix
RADAMPZ	Matrix	N x N equivalent viscous damping ratios
RADAMPG	Matrix	N X N equivalent structural damping ratio

LAMA	LAMA	Eigenvalue summary table
OGPWG	OGPWG	Mass properties output
OQMG1	OQMG	Modal MPC forces
RANCONS	ORGY1	Constraint mode element strain energy table
RANEATC	ORGY1	Attachment mode element strain energy table
RAGCONS	OGPFB	Constraint mode grid point force table
RAGEATC	OGPFB	Attachment mode grid point force table
RAPCONS	OES	Constraint mode ply stress table
RAPEATC	OES	Attachment mode ply stress table
RASCONS	OES	Constraint mode element stress table
RAECONS	OES	Constraint mode element strain table
RASEATC	OES	Attachment mode element stress table
RAEEATC	OES	Attachment mode element strain table
OES1C	OES	Modal Element Stress Table
OES1X	OES	Modal Element Stress Table
OSTR1C	OES	Modal Element Strain Table
OSTR1X	OSTR	Modal Element Strain Table
RAQCONS	OUG	Constraint mode MPC force table
RADCONS	OUG	Constraint mode displacement table
RADEFFM	OUG	Effective inertia displacement table
RAQEATC	OUG	Attachment mode  MPC force table
RADEATC	OUG	Attachment mode displacement table
OUGV1	OUG	Eigenvector Table
RAFCONS	OEF	Constraint mode element force table
RAFEATC	OEF	Attachment mode element force table
OEF1X	OEF	Modal Element Force Table
OGPFB1	OGPFB	Modal Grid Point Force Table
ONRGY1	ONRGY1	Modal Element Strain Energy Table
ONRGY2	ONRGY1

#--------------------

RADCONS - Displacement Constraint Mode
RADDATC - Displacement Distributed Attachment Mode
RADNATC - Displacement Nodal Attachment Mode
RADEATC - Displacement Equivalent Inertia Attachment Mode
RADEFFM - Displacement Effective Inertia Mode

RAECONS - Strain Constraint Mode
RAEDATC - Strain Distributed Attachment Mode
RAENATC - Strain Nodal Attachment Mode
RAEEATC - Strain Equivalent Inertia Attachment Mode

RAFCONS - Element Force Constraint Mode
RAFDATC - Element Force Distributed Attachment Mode
RAFNATC - Element Force Nodal Attachment Mode
RAFEATC - Element Force Equivalent Inertia Attachment Mode

RALDATC - Load Vector Used to Compute the Distributed Attachment M

RANCONS - Strain Energy Constraint Mode
RANDATC - Strain Energy Distributed Attachment Mode
RANNATC - Strain Energy Nodal Attachment Mode
RANEATC - Strain Energy Equivalent Inertia Attachment Mode

RAQCONS - Ply Strains Constraint Mode
RAQDATC - Ply Strains Distributed Attachment Mode
RAQNATC - Ply Strains Nodal Attachment Mode
RAQEATC - Ply Strains Equivalent Inertia Attachment Mode

RARCONS - Reaction Force Constraint Mode
RARDATC - Reaction Force Distributed Attachment Mode
RARNATC - Reaction Force Nodal Attachment Mode
RAREATC - Reaction Force Equivalent Inertia Attachment Mode

RASCONS - Stress Constraint Mode
RASDATC - Stress Distributed Attachment Mode
RASNATC - Stress Nodal Attachment Mode
RASEATC - Stress Equivalent Inertia Attachment Mode

RAPCONS - Ply Stresses Constraint Mode
RAPDATC - Ply Stresses Distributed Attachment Mode
RAPNATC - Ply Stresses Nodal Attachment Mode
RAPEATC - Ply Stresses Equivalent Inertia Attachment Mode

RAGCONS - Grid Point Forces Constraint Mode
RAGDATC - Grid Point Forces Distributed Attachment Mode
RAGNATC - Grid Point Forces Nodal Attachment Mode
RAGEATC - Grid Point Forces Equivalent Inertia Attachment Mode

RADEFMP - Displacement PHA^T * Effective Inertia Mode

RADAMPZ - Viscous Damping Ratio Matrix
RADAMPG - Structural Damping Ratio Matrix

RAFGEN  - Generalized Forces Matrix
BHH     - Modal Viscous Damping Matrix
K4HH    - Modal Structural Damping Matrix
"""
GEOM_TABLES = MSC_GEOM_TABLES + NX_GEOM_TABLES

AUTODESK_MATRIX_TABLES = [
    #b'BELM',
    b'KELM',
    #b'MELM',
] # type: List[bytes]
# this will be split later
TEST_MATRIX_TABLES = [b'ATB', b'BTA', b'MYDOF']

RESULT_TABLES = NX_RESULT_TABLES + MSC_RESULT_TABLES
MATRIX_TABLES = NX_MATRIX_TABLES + MSC_MATRIX_TABLES + AUTODESK_MATRIX_TABLES + TEST_MATRIX_TABLES + [b'MEFF']

#GEOM_TABLES = MSC_GEOM_TABLES
#RESULT_TABLES = MSC_RESULT_TABLES
#MATRIX_TABLES = MSC_MATRIX_TABLES

# TODO: these are weird...
#   RPOSTS1, MAXRATI, RESCOMP, PDRMSG
INT_PARAMS_1 = {
    b'OMODES', b'LGSTRN', b'ADJFRQ', b'BSHDMP', b'BSHDMP4',

    b'POST', b'OPPHIPA', b'OPPHIPB', b'GRDPNT', b'RPOSTS1', b'BAILOUT',
    b'COUPMASS', b'CURV', b'INREL', b'MAXRATI', b'OG',
    b'S1AM', b'S1M', b'DDRMM', b'MAXIT', b'PLTMSG', b'LGDISP', b'NLDISP',
    b'OUNIT2K', b'OUNIT2M', b'RESCOMP', b'PDRMSG', b'LMODES', b'USETPRT',
    b'NOCOMPS', b'OPTEXIT', b'RSOPT', b'GUSTAERO', b'MPTUNIT',
    b'USETSEL', b'NASPRT', b'DESPCH', b'DESPCH1', b'COMPARE', b'DBNBLKS', b'NEWSEQ', b'OLDSEQ',
    b'METHCMRS', b'NOFISR', b'KGGCPCH', b'ERROR', b'DBCDIAG', b'GPECT', b'LSTRN',
    b'DBDROPT', b'SEOP2CV', b'IRES', b'SNORMPRT', b'DBDRNL', b'VMOPT',
    b'OSWPPT', b'KDAMP', b'KDAMPFL', b'MATNL', b'MPCX', b'GEOMPLT', b'NOELOP',
    b'NOGPF', b'PROUT', b'SUPER', b'LGDIS', b'EST', b'SEP1XOVR',
    b'FRSEID', b'HRSEID', b'LRSEID', b'MODACC', b'XFLAG', b'TSTATIC',
    b'NASPDV', b'RMXCRT', b'RMXTRN', b'DBCLEAN', b'LANGLE', b'SEMAPPRT',
    b'FIXEDB', b'AMGOK', b'ASING', b'CNSTRT', b'CURVPLOT', b'CYCIO',
    b'CYCSEQ', b'DBDICT', b'DBINIT', b'DBSET1', b'DBSET2', b'DBSET3', b'DBSET4',
    b'DBSORT', b'DOPT', b'FACTOR', b'ALTSHAPE', b'MODTRK', b'IFTM', b'INRLM',
    b'KINDEX', b'KMIN', b'KMAX', b'LARGEDB', b'LOADINC', b'LOADING', b'LOOP',
    b'LOOPID', b'MODEL', b'MOREK', b'NEWDYN', b'NFECI', b'NINTPTS',
    b'NLAYERS', b'NOELOF', b'NOMSGSTR', b'NONCUP', b'NUMOUT', b'NUMOUT1', b'NUMOUT2',
    b'OPGTKG', b'OPPHIB', b'OUTOPT', b'PKRSP', b'RSPECTRA', b'RSPRINT',
    b'S1G', b'SCRSPEC', b'SEMAPOPT', b'SEQOUT', b'SESEF', b'SKPAMG', b'SKPAMP',
    b'SLOOPID', b'SOLID', b'SPCGEN', b'SRTELTYP', b'SRTOPT', b'START', b'SUBID',
    b'SUBSKP', b'TABID', b'TESTNEG', b'BDMNCON', b'FRUMIN',

    # not defined in qrg...
    b'NT', b'PNCHDB', b'DLOAD', b'NLOAD', b'NOAP', b'NOCMPFLD', b'NODATA',
    b'NODJE', b'NOMECH', b'NOSDR1', b'NOSHADE', b'NOSORT1', b'NOTRED',
    b'NSEGS', b'OLDELM', b'OPADOF', b'OUTPUT', b'P1', b'P2', b'P3', b'PCHRESP',
    b'PLOT', b'PLOTSUP', b'PRTPCH', b'RADLIN', b'RESDUAL', b'S1', b'SDATA',
    b'SEFINAL', b'SEMAP1', b'SKPLOAD', b'SKPMTRX', b'SOLID1', b'SSG3',
    b'PEDGEP', b'ACMSPROC', b'ACMSSEID', b'ACOUS', b'ACOUSTIC', b'ADJFLG',
    b'ADJLDF', b'AEDBCP', b'AESRNDM', b'ARCSIGNS', b'ATVUSE', b'BADMESH', b'BCHNG',
    b'BCTABLE', b'ROTCSV', b'ROTGPF', b'BEARDMP', b'BEARFORC', b'OP2FMT',
    b'LRDISP',

    # ???
    b'CHKSEC', b'CMSMETH', b'CNTNSUB', b'CNTSTPS', b'CONCHG', b'CP',
    b'DBDRPRJ', b'DBDRVER', b'DDAMRUN', b'DESCONX', b'DESEIG', b'DESFINAL',
    b'DESMAX', b'DESSOLAP', b'DIAGOPT',

    b'DOBUCKL', b'DOF123', b'DOMODES', b'DOSTATIC', b'DOTRIP', b'DRESP', b'DSGNOPTX',
    b'DYNAMICX', b'EBULK', b'EIGNFREQ', b'ELOOPID',
    b'FDEPCB', b'FLUIDMP', b'FLUIDSE', b'FMODE', b'FREQDEP', b'FREQDEPS',
    b'GENEL', b'GEOMFLAG', b'GEOMU', b'GKCHNG', b'GLUSET', b'GMCONV', b'GNSTART',
    b'GOODVER', b'GOPH2', b'GRIDFMP', b'GRIDMP', b'HNNLK', b'ICTASET', b'IFPCHNG',
    b'INEP', b'INP2FMT', b'INP4FMT', b'INREL0', b'ITAPE', b'ITOITCNT',
    b'ITOMXITR', b'ITONDVAR', b'ITONGHBR', b'ITONOBJF', b'ITOOPITR', b'ITOPALG',
    b'ITOPALLR', b'ITOPDIAG', b'ITOPOPT', b'ITOSIMP',
    b'IUNIT', b'K4CHNG', b'KCHNG', b'KREDX', b'LANGLES',
    b'LBEARING', b'LDSTI1', b'LMDYN', b'LMODESFL', b'LMSTAT', b'LNUMROT',
    b'LOADGENX', b'LOADREDX', b'LOADU', b'LODCHG', b'LROTOR', b'LTOPOPT',
    b'LUSET', b'LUSETD', b'LUSETS', b'LUSETX', b'MATGENX',
    b'MAXITER', b'MAXSEIDX', b'MBDIFB', b'MBDIFO', b'MBDLMN',
    b'MCHNG', b'MDOF', b'MDTRKFLG', b'MELPG', b'MGRID', b'MLTIMSTR', b'MODESX',
    b'MODETRAK', b'MPIFRHD', b'MPNFLG', b'MREDX', b'MSCOP2', b'NACEXTRA',
    b'NCNOFFST', b'NDISOFP', b'NDVAR', b'NEWSET', b'NGELS', b'NJ', b'NK',
    b'NLBEAR', b'NLCBFOR', b'NMLOOP', b'NMSOL', b'NOA', b'NOASET', b'NOCOMP',
    b'NOFASET', b'NOFGSET', b'NOGENL', b'NOGEOM3', b'NOK4GG', b'NOK4JJ',
    b'NOKGGX', b'NOKJJX', b'NOLSET', b'NOMGG', b'NOMGGX', b'NOMJJX', b'NOQSET',
    b'NORADMAT', b'NORBM', b'NOSE', b'NOSIMP', b'NOSSET', b'NOUE', b'NOUP',
    b'NOYSET', b'NOZSET', b'NQSET', b'NR1OFFST', b'NR2OFFST', b'NR3OFFST',
    b'NROTORS', b'NSE', b'NSKIP0', b'NSOL', b'NSOLF', b'NUMPAN', b'NX',
    b'O2E', b'OADPMAX', b'OALTSHP', b'ODESMAX', b'ODSFLG', b'OMAXR',
    b'OP2SE', b'OP4FMT', b'OP4SE', b'OPGEOM', b'OPTIFCS',
    b'OPTII231', b'OPTII408', b'OPTII411', b'OPTII420', b'OPTIIDMP', b'OPTISNS',
    b'OTAPE', b'OUNIT1', b'OUNIT2', b'OUNIT2R', b'OUTFMP', b'OUTSMP', b'PANELMP',
    b'PBCONT', b'PCHNG', b'PKLLR', b'POSTU', b'PRTMAT', b'PSLGDVX',
    b'PSLOAD', b'PSORT', b'PVALINIT', b'PVALLAST', b'PVALLIST', b'PYCHNG',
    b'REFOPT', b'RESLTOPT', b'RESPSENX',
    b'RMXPANEL', b'ROTPRES', b'ROTPRT', b'RPDFRD', b'RVCHG', b'RVCHG1',
    b'RVCHG2', b'S1AG', b'SAVERSTL', b'SDSRFLAG', b'SEBULK',
    b'SEDMP231', b'SEDMP265', b'SEDMP408', b'SEDMP411', b'SEDMP445', b'SEDMPFLG',
    b'SELDPRS', b'SKIPSE', b'SNDSEIDX', b'SOLFINAL',
    b'SOLNLX', b'SOLNX', b'SOLVSUB', b'SPLINE', b'STOP0', b'STRUCTMP', b'SWEXIST',
    b'TORSIN', b'UACCEL', b'UNIQIDS', b'VUELJUMP', b'VUENEXT',
    b'VUGJUMP', b'VUGNEXT', b'WRTMAT',
    b'XNTIPS', b'XRESLTOP', b'XSEMEDIA', b'XSEUNIT', b'XTIPSCOL',
    b'XYUNIT', b'XZCOLLCT', b'Z2XSING',
    b'ZUZRI1', b'ZUZRI2', b'ZUZRI3', b'ZUZRI4', b'ZUZRI5', b'ZUZRI6', b'ZUZRI7', b'ZUZRI8', b'ZUZRI9', b'ZUZRI10',
    b'ZUZRL1', b'ZUZRL2', b'ZUZRL3', b'ZUZRL4', b'ZUZRL5', b'ZUZRL6', b'ZUZRL7', b'ZUZRL8', b'ZUZRL9', b'ZUZRL10',
    b'DBCPAE', b'DBCPATH',
    b'EXTBEMI', b'EXTBEMO', b'EXTDRUNT', b'EXTUNIT',
    b'UZROLD',
    b'HIRES'
    # no
    #b'SEPS', b'SMALLQ', b'FEPS',
}
FLOAT_PARAMS_1 = {
    b'EPPRT',
    b'WTMASS', b'SNORM', b'PATVER', b'MAXRATIO', b'EPSHT',
    b'SIGMA', b'TABS', b'AUNITS', b'BOLTFACT', b'LMSCAL',
    'DSZERO', b'G', b'GFL', b'LFREQ', b'HFREQ', b'ADPCON',
    b'W3', b'W4', b'W3FL', b'W4FL', b'PREFDB',
    b'EPZERO', b'DSZERO', b'TINY', b'TOLRSC',
    b'FRSPD', b'HRSPD', b'LRSPD', b'MTRFMAX', b'ROTCMRF', b'MTRRMAX',
    b'LAMLIM', b'BIGER', b'BIGER1', b'BIGER2', b'CLOSE',
    b'EPSBIG', b'EPSMALC', b'EPSMALU', b'KDIAG', b'MACH', b'VREF',
    b'STIME', b'TESTSE', b'LFREQFL', b'Q', b'ADPCONS', b'AFNORM', b'AFZERO',
    b'GE', b'MASSDENS',

    # should this be FLOAT_PARAMS_1???
    b'HFREQFL',

    # not defined
    b'CNTSCL',
    b'PRPA', b'PRPHIVZ', b'PRPJ', b'PRRULV', b'RMAX', b'ARF', b'BOV',
    b'ARS', # b'BSHDAMP',
    b'EPSRC',

    # floats - not verified
    b'THRSHOLD', b'SEPS', b'SMALLQ', b'FEPS',
    b'DSNOKD',

    # or integer (not string)
    b'CONFAC',
    b'DFREQ', b'DFRSPCF', b'DSTSPCF', b'DTRSPCF',
    b'DUCTFMAX',
    b'EXTDONE',
    b'FZERO', b'LMFACT', b'MPCZERO',
    b'RESVPGF', b'RESVRAT', b'SWPANGLE', b'UPFAC',
    b'ITODENS',
    b'ITOPCONV',
    b'ITORMAS',
    b'ITOSIMP1',
    b'ITOSIMP2',
    b'MAXRPM',
    b'OBJIN',
    b'PITIME',
    b'RGBEAMA', b'RGBEAME', b'RGLCRIT', b'RGSPRGK',
    b'VOL', b'VOLS',
    b'WGT', b'WGTS',
    b'XSMALLQ',
    b'XUPFAC',
    b'ZUZRR1', b'ZUZRR2', b'ZUZRR3', b'ZUZRR4', b'ZUZRR5', b'ZUZRR6', b'ZUZRR7', b'ZUZRR8', b'ZUZRR9', b'ZUZRR10',
    b'K6ROT',

    # models/msc/units_mass_spring_damper.op2
    b'RBTR',
}
FLOAT_PARAMS_2 = {
    b'BETA', b'CB1', b'CB2', b'CK1', b'CK2', b'CK3', b'CK41', b'CK42',
    b'CM1', b'CM2',
    b'G1', b'G2', b'G3', b'G4', b'G5', b'G6', b'G7', b'G8', b'G9', b'G10',
    b'G11', b'G12', b'G13', b'G14', b'G15', b'G16', b'G17', b'G18', b'G19',
    b'ALPHA1', b'ALPHA2',
    b'CA1', b'CA2',
    b'CP1', b'CP2',
    b'LOADFACS',
    b'ZUZRC1', b'ZUZRC2', b'ZUZRC3', b'ZUZRC4', b'ZUZRC5', b'ZUZRC6', b'ZUZRC7', b'ZUZRC8', b'ZUZRC9', b'ZUZRC10',


    # should this be FLOAT_PARAMS_1???
    #b'EPPRT',
}
INT_PARAMS_2 = {
    b'LOADFACS',
    b'ZUZRC1', b'ZUZRC2', b'ZUZRC3', b'ZUZRC4', b'ZUZRC5', b'ZUZRC6', b'ZUZRC7', b'ZUZRC8', b'ZUZRC9', b'ZUZRC10',
}
#DOUBLE_PARAMS_1 = [] # b'Q'
SATK_STR_PARAMS1 = {
    b'RDRBE',
}
STR_PARAMS_1 = SATK_STR_PARAMS1 | {
    b'POSTEXT', b'PRTMAXIM', b'AUTOSPC', b'OGEOM', b'PRGPST',
    b'RESVEC', b'RESVINER', b'ALTRED', b'OGPS', b'OIBULK', b'OMACHPR',
    b'UNITSYS', b'F56', b'OUGCORD', b'OGEM', b'EXTSEOUT',
    b'CDIF', b'SUPAERO', b'RSCON', b'AUTOMPC', b'DBCCONV',
    b'AUTOSPRT', b'PBRPROP', b'OMID', b'HEATSTAT', b'SECOMB', b'ELEMITER',
    b'ELITASPC', b'DBCONV', b'SHLDAMP', b'COMPMATT', b'SPCSTR', b'ASCOUP',
    b'PRTRESLT', b'SRCOMPS', b'CHECKOUT', b'SEMAP', b'AESMETH', b'RESVALT',
    b'ROTSYNC', b'SYNCDAMP', b'PRGPOST', b'WMODAL', b'SDAMPUP',
    b'COLPHEXA', b'CHKOUT', b'CTYPE', b'DBNAME', b'VUHEXA', b'VUPENTA', b'VUTETRA',
    b'MESH', b'OPTION', b'PRINT', b'SENAME', b'MECHFIX', b'RMXTRAN', b'FLEXINV',
    b'ADSTAT', b'ACOUT', b'ACSYM', b'ACTYPE', b'ADBX', b'AUTOSEEL',
    b'RDSPARSE',
    b'SPARSEDR',
    b'BSHDAMP',
    b'CORROPT',
    b'DBACOUS',
    b'DBALLNOQ',
    b'DBALLX',
    b'DBAPI',
    b'DBAPP',
    b'DBCNT',
    b'DBCOVWRT',
    b'DBDNOPT',
    b'DBDNR', b'DBDNR1', b'DBDNX', b'DBEXT', b'DBGOA', b'DBMAP',
    b'DBOFP2X', b'DBOFPX', b'DBRCVX', b'DBSCRR', b'DBUPOPT', b'DBUPR',
    b'DBUPX', b'DBXSEDR', b'DBXSEDRR', b'DBZUZR', b'DSOR', b'DSOX',
    b'DVGRDN', b'DYNSPCF', b'EQVSCR', b'EXTDROUT',
    b'FLEXINCR', b'FTL', b'GDAMPF', b'GEOCENT', b'IFPSCR', b'IFPSOPT',
    b'IFPX', b'IFPXOPT', b'MASTER', b'MODEOUT',
    b'NXVER', b'OAPP', b'OCMP', b'OEE', b'OEEX', b'OEF', b'OEFX', b'OEPT',
    b'OES', b'OESE', b'OESX', b'OGPF', b'OMPT', b'OPG', b'OPTIM', b'OQG',
    b'OUG', b'OUMU', b'OUTSCR', b'PANAME', b'QSETREM', b'RESVSE', b'RESVSLI',
    b'RESVSO', b'RSATT', b'SAVEOFP', b'SAVERST', b'SCRATCH', b'SDRPOPT',
    b'SECOMB0', b'SELRNG', b'SERST', b'SOFTEXIT', b'SOLAPPI', b'SOLTYPI',
    b'TDB0', b'TDBX', b'UPDTBSH',
    b'USETSTR1', b'USETSTR2', b'USETSTR3', b'USETSTR4',
    b'VMOPTSET', b'VUBEAM', b'VUQUAD4', b'VUTRIA3', b'WRN', b'XAUTOSPT',
    b'XRESVECA', b'XRESVECO', b'XRESVIRA', b'XRESVIRO',
    b'ZUZRCL1', b'ZUZRCL2', b'ZUZRCL3', b'ZUZRCL4', b'ZUZRCL5', b'ZUZRCL6', b'ZUZRCL7', b'ZUZRCL8', b'ZUZRCL9', b'ZUZRCL10',
    b'ZUZRCH1', b'ZUZRCH2', b'ZUZRCH3', b'ZUZRCH4', b'ZUZRCH5', b'ZUZRCH6', b'ZUZRCH7', b'ZUZRCH8', b'ZUZRCH9', b'ZUZRCH10',
    b'APPI', b'APPF',

    # part of param, checkout
    b'PRTBGPDT', b'PRTCSTM', b'PRTEQXIN', b'PRTGPDT',
    b'PRTGPL', b'PRTGPTT', b'PRTMGG', b'PRTPG',

    # superelements
    b'EXTOUT', b'SESDAMP',

    # TODO: remove these as they're in the matrix test and are user
    #       defined PARAMs; arguably all official examples should just work
    # TODO: add an option for custom PARAMs
    b'ADB', b'AEDB', b'MREDUC', b'OUTDRM', b'OUTFORM', b'REDMETH', b'DEBUG',
    b'AEDBX', b'AERO', b'AUTOSUP0', b'AXIOPT',
}
def _check_unique_sets(*sets: List[set[str]]):
    """verifies that the sets are unique"""
    for i, seti in enumerate(sets):
        for unused_j, setj in enumerate(sets[i+1:]):
            intersectioni = seti.intersection(setj)
            assert len(intersectioni) == 0, intersectioni

_check_unique_sets(INT_PARAMS_1, FLOAT_PARAMS_1, FLOAT_PARAMS_2, STR_PARAMS_1)


class OP2_Scalar(OP2Common, FortranFormat):
    """Defines an interface for the Nastran OP2 file."""
    @property
    def total_effective_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['EFMFSMS']

    @property
    def effective_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['EFMASSS']

    @property
    def rigid_body_mass_matrix(self):
        """6x6 matrix"""
        return self.matrices['RBMASS']

    @property
    def modal_effective_mass_fraction(self):
        """6xnmodes matrix"""
        return self.matrices['EFMFACS']#.dataframe

    @property
    def modal_participation_factors(self):
        """6xnmodes matrix"""
        return self.matrices['MPFACS']#.dataframe

    @property
    def modal_effective_mass(self):
        """6xnmodes matrix"""
        return self.matrices['MEFMASS']#.dataframe

    @property
    def modal_effective_weight(self):
        """6xnmodes matrix"""
        return self.matrices['MEFWTS']#.dataframe


    @property
    def monitor1(self):
        self.deprecated('op2.monitor1', 'op2.op2_results.monitor1', '1.4')
        return self.op2_results.monitor1
    @monitor1.setter
    def monitor1(self, monitor1):
        self.deprecated('op2.monitor1', 'op2.op2_results.monitor1', '1.4')
        self.op2_results.monitor1 = monitor1

    @property
    def monitor3(self):
        self.deprecated('op2.monitor3', 'op2.op2_results.monitor3', '1.4')
        return self.op2_results.monitor3
    @monitor3.setter
    def monitor3(self, monitor3):
        self.deprecated('op2.monitor3', 'op2.op2_results.monitor3', '1.4')
        self.op2_results.monitor3 = monitor3

    @property
    def matrix_tables(self):
        return MATRIX_TABLES

    def set_as_nx(self):
        self.is_nx = True
        self.is_msc = False
        self.is_autodesk = False
        self.is_nasa95 = False
        self.is_optistruct = False
        self._nastran_format = 'nx'

    def set_as_msc(self):
        self.is_nx = False
        self.is_msc = True
        self.is_autodesk = False
        self.is_nasa95 = False
        self.is_optistruct = False
        self._nastran_format = 'msc'

    def set_as_autodesk(self):
        self.is_nx = False
        self.is_msc = False
        self.is_autodesk = True
        self.is_nasa95 = False
        self.is_optistruct = False
        self._nastran_format = 'autodesk'

    def set_as_nasa95(self):
        self.is_nx = False
        self.is_msc = False
        self.is_autodesk = False
        self.is_optistruct = False
        self.is_nasa95 = True
        self._nastran_format = 'nasa95'
        self.reader_oes._read_oes1_loads = self.reader_oes._read_oes1_loads_nasa95
        self.reader_oef._read_oef1_loads = self.reader_oef._read_oef1_loads_nasa95

        if hasattr(self, 'reader_geom2') and hasattr(self.reader_geom2, '_read_cquad4_nasa95'):
            self.reader_geom2.geom2_map[(5408, 54, 261)] = ['CQUAD4', self.reader_geom2._read_cquad4_nasa95]

    def set_as_optistruct(self):
        self.is_nx = False
        self.is_msc = False
        self.is_autodesk = False
        self.is_nasa95 = False
        self.is_optistruct = True
        self._nastran_format = 'optistruct'

    def __init__(self, debug=False, log=None, debug_file=None):
        """
        Initializes the OP2_Scalar object

        Parameters
        ----------
        debug : bool/None; default=True
            used to set the logger if no logger is passed in
                True:  logs debug/info/warning/error messages
                False: logs info/warning/error messages
                None:  logs warning/error messages
        log : Log()
            a logging object to write debug messages to
            (.. seealso:: import logging)
        debug_file : str; default=None (No debug)
            sets the filename that will be written to

        """
        assert debug is None or isinstance(debug, bool), 'debug=%r' % debug

        self.log = get_logger2(log, debug=debug, encoding='utf-8')
        self._count = 0
        self.op2_filename = None
        self.bdf_filename = None
        self.f06_filename = None
        self.des_filename = None
        self.h5_filename = None
        self._encoding = 'utf8'

        #: should a MATPOOL "symmetric" matrix be stored as symmetric
        #: it takes double the RAM, but is easier to use
        self.apply_symmetry = True
        OP2Common.__init__(self)

        FortranFormat.__init__(self)

        self.is_vectorized = False
        self._close_op2 = True

        self.result_names = set()

        self.grid_point_weight = {}
        self.words = []
        self.debug = debug
        self._last_comment = None
        #self.debug = True
        #self.debug = False
        #debug_file = None
        if debug_file is None:
            self.debug_file = None
        else:
            assert isinstance(debug_file, str), debug_file
            self.debug_file = debug_file

        self.op2_reader = OP2Reader(self)

    def set_subcases(self, subcases=None):
        """
        Allows you to read only the subcases in the list of isubcases

        Parameters
        ----------
        subcases : List[int, ...] / int; default=None->all subcases
            list of [subcase1_ID,subcase2_ID]

        """
        #: stores the set of all subcases that are in the OP2
        #self.subcases = set()
        if subcases is None or subcases == []:
            #: stores if the user entered [] for isubcases
            self.is_all_subcases = True
            self.valid_subcases = []
        else:
            #: should all the subcases be read (default=True)
            self.is_all_subcases = False

            if isinstance(subcases, int):
                subcases = [subcases]

            #: the set of valid subcases -> set([1,2,3])
            self.valid_subcases = set(subcases)
        self.log.debug(f'set_subcases - subcases = {self.valid_subcases}')

    def set_transient_times(self, times):  # TODO this name sucks...
        """
        Takes a dictionary of list of times in a transient case and
        gets the output closest to those times.

        Examples
        --------
        >>> times = {subcase_id_1: [time1, time2],
                     subcase_id_2: [time3, time4]}

        .. warning:: I'm not sure this still works...

        """
        expected_times = {}
        for (isubcase, etimes) in times.items():
            etimes = list(times)
            etimes.sort()
            expected_times[isubcase] = array(etimes)
        self.expected_times = expected_times

    def _get_table_mapper(self):
        """gets the dictionary of function3 / function4"""

        # MSC table mapper
        reader_lama = self.reader_lama
        reader_ogpwg = self.reader_ogpwg
        reader_ogpf = self.reader_ogpf
        reader_onr = self.reader_onr
        reader_opg = self.reader_opg
        reader_opr = self.reader_opr
        reader_oqg = self.reader_oqg
        reader_ogs = self.reader_ogs

        # oug
        reader_oug = self.reader_oug
        reader_ougpk = self.reader_ougpk
        reader_otemp = self.reader_otemp

        # oef
        reader_oef = self.reader_oef
        reader_oefpk = self.reader_oefpk

        # oes
        reader_oes = self.reader_oes

        satk_tables = {
            b'OUGPK1'  : [self.reader_ougpk._read_ougpk1_3, self.reader_ougpk._read_ougpk1_4],
            b'OEFPK1'  : [self.reader_oefpk._read_oefpk1_3, self.reader_oefpk._read_oefpk1_4],
        }

        nasa95_tables = {
        b'OESC1'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
        }

        table_mapper_geometry = {
            # -----------------------------------------------------------
            # geometry
            b'GEOM1' : [self._table_passer, self._table_passer], # GEOM1-Geometry-related bulk data
            b'GEOM2' : [self._table_passer, self._table_passer], # GEOM2-element connectivity and SPOINT-related data
            b'GEOM3' : [self._table_passer, self._table_passer], # GEOM3-Static and thermal loads
            b'GEOM4' : [self._table_passer, self._table_passer], # GEOM4-constraints, DOF membership entries, MPC, and R-type element data

            # superelements
            b'GEOM1S' : [self._table_passer, self._table_passer],  # GEOMx + superelement
            b'GEOM2S' : [self._table_passer, self._table_passer],
            b'GEOM3S' : [self._table_passer, self._table_passer],
            b'GEOM4S' : [self._table_passer, self._table_passer],

            b'GEOM1VU' : [self._table_passer, self._table_passer],
            b'GEOM2VU' : [self._table_passer, self._table_passer],
            b'BGPDTVU' : [self._table_passer, self._table_passer],

            b'GEOM1N' : [self._table_passer, self._table_passer],
            b'GEOM2N' : [self._table_passer, self._table_passer],
            b'GEOM3N' : [self._table_passer, self._table_passer],
            b'GEOM4N' : [self._table_passer, self._table_passer],

            b'GEOM1OLD' : [self._table_passer, self._table_passer],
            b'GEOM2OLD' : [self._table_passer, self._table_passer],
            b'GEOM3OLD' : [self._table_passer, self._table_passer],
            b'GEOM4OLD' : [self._table_passer, self._table_passer],

            b'EPT' : [self._table_passer, self._table_passer],  # elements
            b'EPTS' : [self._table_passer, self._table_passer],  # elements - superelements
            b'EPTOLD' : [self._table_passer, self._table_passer],

            b'MPT' : [self._table_passer, self._table_passer],  # materials
            b'MPTS' : [self._table_passer, self._table_passer],  # materials - superelements

            b'DYNAMIC' : [self._table_passer, self._table_passer],
            b'DYNAMICS' : [self._table_passer, self._table_passer],
            b'DIT' : [self._table_passer, self._table_passer],
            b'DITS' : [self._table_passer, self._table_passer],

            # this comment may refer to CSTM?
            #F:\work\pyNastran\examples\Dropbox\pyNastran\bdf\cards\test\test_mass_01.op2
            #F:\work\pyNastran\examples\matpool\gpsc1.op2
            b'AXIC': [self._table_passer, self._table_passer],

            # EDT - aero cards
            #       element deformation, aerodynamics, p-element, divergence analysis,
            #       and iterative solver input (includes SET1 entries)
            b'EDT' : [self._table_passer, self._table_passer],
            b'EDTS' : [self._table_passer, self._table_passer],

            # contact/glue
            b'CONTACT' : [self._table_passer, self._table_passer],
            b'CONTACTS' : [self._table_passer, self._table_passer],

            b'EDOM' : [self._table_passer, self._table_passer],     # optimization
            b'VIEWTB' : [self._table_passer, self._table_passer],   # view elements
            # =========================end geom passers=========================
        }

        table_mapper = {
            # per NX
            b'OESVM1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],    # isat_random
            b'OESVM1C' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],   # isat_random
            b'OSTRVM1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],   # isat_random
            b'OSTRVM1C' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # isat_random

            b'OSTRVM2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],

            b'OESVM2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],    # big random
            b'OES2C' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            b'OSTR2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4], # TODO: disable
            b'OSTR2C' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],

            b'OTEMP1' : [reader_otemp._read_otemp1_3, reader_otemp._read_otemp1_4],
            # --------------------------------------------------------------------------
            # MSC TABLES
            # common tables

            # unorganized
            b'RADCONS': [reader_oug._read_oug1_3, reader_oug._read_oug_4], # Displacement Constraint Mode (OUG)
            b'RADEFFM': [reader_oug._read_oug1_3, reader_oug._read_oug_4], # Displacement Effective Inertia Mode (OUG)
            b'RADEATC': [reader_oug._read_oug1_3, reader_oug._read_oug_4], # Displacement Equivalent Inertia Attachment mode (OUG)

            # broken - isat_launch_100hz.op2 - wrong numwide
            # spc forces
            b'RAQCONS': [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4], # Constraint mode MPC force table (OQG)
            b'RAQEATC': [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4], # Attachment mode MPC force table (OQG)
            #b'RAQCONS': [self._table_passer, self._table_passer], # temporary
            #b'RAQEATC': [self._table_passer, self._table_passer], # temporary

            b'RAREATC': [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4], # spc forces?
            b'RARCONS': [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4], # spc forces?

            # element forces
            b'RAFCONS': [reader_oef._read_oef1_3, reader_oef._read_oef1_4], # Element Force Constraint Mode (OEF)
            b'RAFEATC': [reader_oef._read_oef1_3, reader_oef._read_oef1_4], # Element Force Equivalent Inertia Attachment mode (OEF)
            #b'RAFCONS': [self._table_passer, self._table_passer], # temporary
            #b'RAFEATC': [self._table_passer, self._table_passer], # temporary

            # grid point forces
            b'RAGCONS': [reader_ogpf._read_ogpf1_3, reader_ogpf._read_ogpf1_4], # Grid Point Forces Constraint Mode (OGPFB)
            b'RAGEATC': [reader_ogpf._read_ogpf1_3, reader_ogpf._read_ogpf1_4], # Grid Point Forces Equivalent Inertia Attachment mode (OEF)
            #b'RAGCONS': [self._table_passer, self._table_passer], # Grid Point Forces Constraint Mode (OGPFB)
            #b'RAGEATC': [self._table_passer, self._table_passer], # Grid Point Forces Equivalent Inertia Attachment mode (OEF)

            # stress
            b'RAPCONS': [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # Constraint mode ply stress table (OES)
            b'RAPEATC': [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # Attachment mode ply stress table (OES)
            #b'RAPCONS': [self._table_passer, self._table_passer], # Constraint mode ply stress table (OES)
            #b'RAPEATC': [self._table_passer, self._table_passer], # Attachment mode ply stress table (OES)

            # stress
            b'RASCONS': [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # Stress Constraint Mode (OES)
            b'RASEATC': [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # Stress Equivalent Inertia Attachment mode (OES)
            #b'RASCONS': [self._table_passer, self._table_passer], # temporary
            #b'RASEATC': [self._table_passer, self._table_passer], # temporary

            # strain
            b'RAEEATC': [reader_oes._read_oes1_3, reader_oes._read_ostr1_4], # Strain Equivalent Inertia Attachment mode (OES)
            b'RAECONS': [reader_oes._read_oes1_3, reader_oes._read_ostr1_4], # Strain Constraint Mode (OSTR)
            #b'RAEEATC': [self._table_passer, self._table_passer], # temporary
            #b'RAECONS': [self._table_passer, self._table_passer], # temporary

            # strain energy
            b'RANEATC' : [reader_onr._read_onr1_3, reader_onr._read_onr1_4], # Strain Energy Equivalent Inertia Attachment mode (ORGY1)
            b'RANCONS': [reader_onr._read_onr1_3, reader_onr._read_onr1_4], # Constraint mode element strain energy table (ORGY1)
            #b'RANEATC': [self._table_passer, self._table_passer], # Strain Energy Equivalent Inertia Attachment mode (ORGY1)
            #b'RANCONS': [self._table_passer, self._table_passer], # Constraint mode element strain energy table (ORGY1)


            #b'TOL': [self._table_passer, self._table_passer],

            b'MATPOOL': [self._table_passer, self._table_passer], # DMIG bulk data entries

            b'RSOUGV1': [self._table_passer, self._table_passer],
            b'RESOES1': [self._table_passer, self._table_passer],
            b'RESEF1' : [self._table_passer, self._table_passer],
            b'DESCYC' : [self._table_passer, self._table_passer],
            #b'AEMONPT' : [self._read_aemonpt_3, self._read_aemonpt_4],
            #=======================
            # OEF
            # element forces
            #b'OEFITSTN' : [self._table_passer, self._table_passer], # works
            b'OEFITSTN' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],
            b'OEFIT' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],  # failure indices
            b'OEF1X' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],  # element forces at intermediate stations
            b'OEF1'  : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],  # element forces or heat flux
            b'HOEF1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],  # element heat flux

            # off force
            b'OEF2' : [reader_oef._read_oef2_3, reader_oef._read_oef2_4],  # element forces or heat flux
            #=======================
            # OQG
            # spc forces
            # OQG1/OQGV1 - spc forces in the nodal frame
            # OQP1 - scaled response spectra - spc-forces
            b'OQG1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQG2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],

            b'OQGV1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGV2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],

            b'PSDF' : [self._read_psdf_3, self._read_psdf_4],  # MSC NASA/goesr

            #=======================
            # OPG
            # applied loads
            b'OPG1'  : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],  # applied loads in the nodal frame
            b'OPGV1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],  # solution set applied loads?
            b'OPNL1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],  # nonlinear loads

            b'OPG2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],   # applied loads in the nodal frame
            b'OPNL2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],  # nonlinear loads

            #=======================
            # OGPFB1
            # grid point forces
            b'OGPFB1' : [reader_ogpf._read_ogpf1_3, reader_ogpf._read_ogpf1_4],  # grid point forces
            #b'OGPFB2' : [reader_ogpf._read_ogpf1_3, reader_ogpf._read_ogpf1_4],  # grid point forces

            #=======================
            # ONR/OEE
            # strain energy density
            b'ONRGY'  : [reader_onr._read_onr1_3, reader_onr._read_onr1_4],
            b'ONRGY1' : [reader_onr._read_onr1_3, reader_onr._read_onr1_4],  # strain energy density
            b'ONRGY2':  [reader_onr._read_onr2_3, reader_onr._read_onr1_4],
            #b'ONRGY2':  [self._table_passer, self._table_passer],
            #===========================================================
            # OES
            # stress
            # OES1C - Table of composite element stresses or strains in SORT1 format
            # OESRT - Table of composite element ply strength ratio. Output by SDRCOMP
            b'OES1X1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # stress - nonlinear elements
            b'OES1'   : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # stress - linear only
            b'OES1X'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # element stresses at intermediate stations & nonlinear stresses
            b'OES1C'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # stress - composite
            b'OESCP'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # stress - nonlinear???
            b'OESRT'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4], # ply strength ratio

            # strain
            b'OSTR1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # strain - autodesk/9zk6b5uuo.op2
            b'OSTR1X'  : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # strain - isotropic
            b'OSTR1C'  : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # strain - composite
            b'OESTRCP' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # strain - composite ply

            b'OSTR1PL' : [self._table_passer, self._table_passer],  # Table of ply strains-plastic in SORT1 format
            b'OSTR1THC' : [self._table_passer, self._table_passer], # ????
            b'OSTR1CR' : [self._table_passer, self._table_passer],  # ????
            #b'OEFIIP'

            # special nonlinear tables
            # OESNLBR - Slideline stresses
            # OESNLXD - Nonlinear transient stresses
            # OESNLXR - Nonlinear stress
            #           Table of nonlinear element stresses in SORT1 format and appended for all subcases

            b'OESNLXR' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # nonlinear stresses
            b'OESNLXD' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # nonlinear transient stresses
            b'OESNLBR' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESNL1X' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],

            b'OESNL2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            b'OESNLXR2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            b'OESNLBR2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            #b'OESNLXR2' : [self._table_passer, self._table_passer],
            #b'OESNLBR2' : [self._table_passer, self._table_passer],

            # off stress
            b'OES2'    : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],  # stress - linear only - disabled; need better tests

            #=======================
            # OUG
            # displacement/velocity/acceleration/eigenvector/temperature
            b'OUG1'    : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # displacements in nodal frame
            b'OVG1'    : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # velocity in nodal frame
            b'OAG1'    : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # accelerations in nodal frame

            b'OUG1F'    : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # acoustic displacements in ? frame

            #b'OVG'   : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # aeroelastic velocity
            b'OUGV1'   : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # displacements in nodal frame
            b'BOUGV1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # OUG1 on the boundary???
            b'BOUGF1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # OUG1 on the boundary???
            b'OUGV1PAT': [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # OUG1 + coord ID
            b'TOUGV1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # grid point temperature
            b'ROUGV1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # relative OUG
            b'OPHSA'   : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # Displacement table in SORT1
            b'OUXY1'   : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # Displacements in SORT1 format for h-set or d-set.
            b'OUGPC1'  : [reader_oug._read_ougpc1_3, reader_oug._read_ougpc_4],  # panel contributions
            b'OUGPC2'  : [reader_oug._read_ougpc2_3, reader_oug._read_ougpc_4],  # panel contributions
            b'OUGF1' : [reader_oug._read_oug1_3, reader_oug._read_oug_4], # Acoustic pressures at microphone points in SORT1 format
            b'OUGF2' : [reader_oug._read_oug2_3, reader_oug._read_oug_4], # Acoustic pressures at microphone points in SORT1 format

            b'OUGV2'   : [reader_oug._read_oug2_3, reader_oug._read_oug_4],  # displacements in nodal frame
            b'ROUGV2'  : [reader_oug._read_oug2_3, reader_oug._read_oug_4],  # relative OUG
            b'OUXY2'   : [reader_oug._read_oug2_3, reader_oug._read_oug_4],  # Displacements in SORT2 format for h-set or d-set.

            #  scaled response spectra - ABS / NRL / SRSS
            b'OUPV1' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],    # displacement, velocity, acceleration
            b'OQP1'  : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],    # spc forces
            b'OQP2'  : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],    # spc forces
            b'DOEF1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],   # forces
            b'DOES1'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # stress
            b'DOSTR1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4], # strain

            # modal contribution
            b'OUGMC1'  : [reader_oug._read_oug1_3, reader_oug._read_ougmc_4],  # displacement
            b'OQGMC1'  : [reader_oqg._read_oqg1_3, reader_oug._read_ougmc_4],  # spc/mpc? force
            b'OESMC1'  : [reader_oes._read_oes1_3, reader_oes._read_oesmc_4],  # stress
            b'OSTRMC1'  : [reader_oes._read_oes1_3, reader_oes._read_oesmc_4], # strain

            #F:\work\pyNastran\examples\Dropbox\move_tpl\sbuckl2a.op2
            # post-buckling / nonlinear buckling
            b'OCRUG' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],   # displacement
            b'OCRPG' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],  # loads

            # eigenvectors
            b'OPHIG' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],    # basic coordinate system
            b'BOPHIG' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],   # basic coordinate system
            b'BOPHIGF' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # basic coordinate system for the fluid portion of the model
            b'BOPHIGS' : [reader_oug._read_oug1_3, reader_oug._read_oug_4],  # basic coordinate system for the structural portion of the model

            b'BOPG1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],  # applied loads in basic coordinate system

            #=======================
            # extreme values of the respective table
            b'OUGV1MX' : [self._table_passer, self._table_passer],
            b'OEF1MX' : [self._table_passer, self._table_passer],
            b'OES1MX' : [self._table_passer, self._table_passer],

            #=======================
            # contact
            b'OQGCF1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4], # Contact force at grid point.
            b'OQGCF2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4], # Contact force at grid point.

            b'OSPDS1' : [reader_oqg._read_opsds1_3, reader_oqg._read_opsds1_4],  # Final separation distance.
            b'OSPDS2' : [self._nx_table_passer, self._table_passer],

            b'OSPDSI1' : [reader_oqg._read_opsdi1_3, reader_oqg._read_opsdi1_4], # Initial separation distance.
            b'OSPDSI2' : [self._nx_table_passer, self._table_passer], # Output contact separation distance results.

            #b'OBC1' : [self._read_obc1_3, self._read_obc1_4],
            #b'OBC2' : [self._nx_table_passer, self._table_passer], # Contact pressures and tractions at grid points.

            #b'OSLIDE1'
            b'OCPSDFC' : [self._nx_table_passer, self._table_passer],
            b'OCCORFC' : [self._nx_table_passer, self._table_passer],

            # Glue normal and tangential tractions at grid point in basic coordinate system
            b'OBG1' : [self._nx_table_passer, self._table_passer],
            b'OBG2' : [self._nx_table_passer, self._table_passer],

            # Glue forces at grid point in basic coordinate system
            b'OQGGF1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGGF2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],

            # Table of Euler Angles for transformation from material to basic coordinate system in the:
            b'TRMBU' : [self._nx_table_passer, self._table_passer], # undeformed configuration
            b'TRMBD' : [self._nx_table_passer, self._table_passer],
            #=======================
            # OGPWG
            # grid point weight
            b'OGPWG'  : [reader_ogpwg._read_ogpwg_3, reader_ogpwg._read_ogpwg_4],  # grid point weight
            b'OGPWGM' : [reader_ogpwg._read_ogpwg_3, reader_ogpwg._read_ogpwg_4],  # modal? grid point weight

            #=======================
            # OGS
            # grid point stresses
            b'OGS1' : [reader_ogs._read_ogs1_3, reader_ogs._read_ogs1_4],  # grid point stresses
            #b'OGS2' : [self._read_ogs1_3, self._read_ogs1_4],  # grid point stresses

            b'OGSTR1' : [reader_ogs._read_ogstr1_3, reader_ogs._read_ogstr1_4],  # grid point strains
            #=======================
            # eigenvalues
            b'BLAMA' : [reader_lama._read_buckling_eigenvalue_3, reader_lama._read_buckling_eigenvalue_4], # buckling eigenvalues
            b'CLAMA' : [reader_lama._read_complex_eigenvalue_3, reader_lama._read_complex_eigenvalue_4],   # complex eigenvalues
            b'LAMA'  : [reader_lama._read_real_eigenvalue_3, reader_lama._read_real_eigenvalue_4],         # eigenvalues
            b'LAMAS' : [reader_lama._read_real_eigenvalue_3, reader_lama._read_real_eigenvalue_4],         # eigenvalues-structure
            b'LAMAF' : [reader_lama._read_real_eigenvalue_3, reader_lama._read_real_eigenvalue_4],         # eigenvalues-fluid

            # ===passers===
            #b'EQEXIN': [self._table_passer, self._table_passer],
            #b'EQEXINS': [self._table_passer, self._table_passer],

            b'GPDT' : [self._table_passer, self._table_passer],     # grid points?
            b'BGPDT' : [self._table_passer, self._table_passer],    # basic grid point definition table
            b'BGPDTS' : [self._table_passer, self._table_passer],
            b'BGPDTOLD' : [self._table_passer, self._table_passer],

            b'PVT' : [self._read_pvto_3, self._read_pvto_4], # PVT - Parameter Variable Table
            b'PVTS' : [self._read_pvto_3, self._read_pvto_4], # ???
            b'PVT0' : [self._read_pvto_3, self._read_pvto_4],  # user parameter value table
            b'TOLD' : [self._table_passer, self._table_passer],
            #b'CASECC' : [self._table_passer, self._table_passer],  # case control deck
            #b'XCASECC' : [self._table_passer, self._table_passer], # ???

            b'STDISP' : [self._table_passer, self._table_passer], # matrix?
            b'AEDISP' : [self._table_passer, self._table_passer], # matrix?
            #b'TOLB2' : [self._table_passer, self._table_passer], # matrix?

            b'FOL' : [self._table_passer, self._table_passer],
            b'PERF' : [self._table_passer, self._table_passer],

            # DSCMCOL - Correlation table for normalized design sensitivity coefficient matrix.
            #           Output by DSTAP2.
            # DBCOPT - Design optimization history table for ???

            b'OEKE1' : [self._table_passer, self._table_passer],
            #b'DSCMCOL' : [self._table_passer, self._table_passer],
            #b'DBCOPT' : [self._table_passer, self._table_passer],
            #b'FRL0': [self._table_passer, self._table_passer],  # frequency response list

            #==================================
            # modal participation factors
            b'OFMPF2M' : [self._read_mpf_3, self._read_mpf_4], # OFMPF2M Table of fluid mode participation factors by normal mode.
            b'OLMPF2M' : [self._read_mpf_3, self._read_mpf_4], # OLMPF2M Load mode participation factors by normal mode.
            b'OPMPF2M' : [self._read_mpf_3, self._read_mpf_4], # OPMPF2M Panel mode participation factors by normal mode.
            b'OSMPF2M' : [self._read_mpf_3, self._read_mpf_4],

            # OGMPF2M Grid mode participation factors by normal mode.
            b'OGPMPF2M' : [self._read_mpf_3, self._read_mpf_4],

            #OFMPF2E Table of fluid mode participation factors by excitation frequencies.
            #OSMPF2E Table of structure mode participation factors by excitation frequencies.
            #OPMPF2E Table of panel mode participation factors by excitation frequencies.
            #OLMPF2E Table of load mode participation factors by excitation frequencies.
            #OGMPF2E Table of grid mode participation factors by excitation frequencies.

            #==================================
            #b'GPL': [self._table_passer, self._table_passer],
            #b'OMM2' : [self._table_passer, self._table_passer],  # max/min table - kinda useless
            b'ERRORN' : [self._table_passer, self._table_passer],  # p-element error summary table
            #==================================
            b'OUG2T' : [self._table_passer, self._table_passer],

            # nx cohesive zone
            b'ODAMGCZT' : [self._nx_table_passer, self._table_passer], # nx - Table of tractions for cohesive elements
            b'ODAMGCZR' : [self._nx_table_passer, self._table_passer], # nx - Table of relative displacements for cohesive elements
            b'ODAMGCZD' : [self._nx_table_passer, self._table_passer], # nx - Table of damage values for cohesive elements

            # Normalized Mass Density
            b'ONMD' : [self.reader_onmd._read_onmd_3, self.reader_onmd._read_onmd_4],
        }
        table_mapper_random = {
            # random OUG (displacement, velocity, acceleration)
            b'OUGATO1' : [reader_oug._read_oug1_3, reader_oug._read_oug_ato],
            b'OUGCRM1' : [reader_oug._read_oug1_3, reader_oug._read_oug_crm],
            b'OUGPSD1' : [reader_oug._read_oug1_3, reader_oug._read_oug_psd],
            b'OUGRMS1' : [reader_oug._read_oug1_3, reader_oug._read_oug_rms],
            b'OUGNO1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_no],

            b'OUGATO2' : [reader_oug._read_oug2_3, reader_oug._read_oug_ato],
            b'OUGCRM2' : [reader_oug._read_oug2_3, reader_oug._read_oug_crm],
            b'OUGPSD2' : [reader_oug._read_oug2_3, reader_oug._read_oug_psd],
            b'OUGRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
            b'OUGNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random
            #b'OUGRMS2' : [reader_oug._read_oug2_3, reader_oug._read_oug_rms],  # buggy on isat random
            #b'OUGNO2'  : [reader_oug._read_oug2_3, reader_oug._read_oug_no],  # buggy on isat random

            # velocity
            b'OVGATO1' : [reader_oug._read_oug1_3, reader_oug._read_oug_ato],
            b'OVGCRM1' : [reader_oug._read_oug1_3, reader_oug._read_oug_crm],
            b'OVGPSD1' : [reader_oug._read_oug1_3, reader_oug._read_oug_psd],
            b'OVGRMS1' : [reader_oug._read_oug1_3, reader_oug._read_oug_rms],
            b'OVGNO1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_no],

            b'OVGATO2' : [reader_oug._read_oug2_3, reader_oug._read_oug_ato],
            b'OVGCRM2' : [reader_oug._read_oug2_3, reader_oug._read_oug_crm],
            b'OVGPSD2' : [reader_oug._read_oug2_3, reader_oug._read_oug_psd],
            #b'OVGRMS2' : [self._table_passer, self._table_passer],
            #b'OVGNO2'  : [self._table_passer, self._table_passer],
            b'OVGRMS2' : [reader_oug._read_oug2_3, reader_oug._read_oug_rms],
            b'OVGNO2'  : [reader_oug._read_oug2_3, reader_oug._read_oug_no],

            # acceleration
            b'OAGATO1' : [reader_oug._read_oug1_3, reader_oug._read_oug_ato],
            b'OAGCRM1' : [reader_oug._read_oug1_3, reader_oug._read_oug_crm],
            b'OAGPSD1' : [reader_oug._read_oug1_3, reader_oug._read_oug_psd],
            b'OAGRMS1' : [reader_oug._read_oug1_3, reader_oug._read_oug_rms],
            b'OAGNO1'  : [reader_oug._read_oug1_3, reader_oug._read_oug_no],

            b'OAGATO2' : [reader_oug._read_oug2_3, reader_oug._read_oug_ato],
            b'OAGCRM2' : [reader_oug._read_oug2_3, reader_oug._read_oug_crm],
            b'OAGPSD2' : [reader_oug._read_oug2_3, reader_oug._read_oug_psd],
            #b'OAGRMS2' : [self._table_passer, self._table_passer],
            #b'OAGNO2'  : [self._table_passer, self._table_passer],
            b'OAGRMS2' : [reader_oug._read_oug2_3, reader_oug._read_oug_rms],
            b'OAGNO2'  : [reader_oug._read_oug2_3, reader_oug._read_oug_no],

            # applied load
            b'OPGATO1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],
            b'OPGCRM1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],
            b'OPGPSD1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],
            b'OPGRMS1' : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],
            b'OPGNO1'  : [reader_opg._read_opg1_3, reader_opg._read_opg1_4],

            b'OPGATO2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],
            b'OPGCRM2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],
            b'OPGPSD2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],
            #b'OPGRMS2' : [self._table_passer, self._table_passer],
            #b'OPGNO2'  : [self._table_passer, self._table_passer],
            b'OPGRMS2' : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],
            b'OPGNO2'  : [reader_opg._read_opg2_3, reader_opg._read_opg1_4],

            # SPC/MPC tables depending on table_code
            # SPC - NX/MSC
            # MPC - MSC
            b'OQGATO1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGCRM1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGPSD1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGRMS1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],
            b'OQGNO1'  : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_4],

            b'OQGATO2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],
            b'OQGCRM2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],
            b'OQGPSD2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],
            b'OQGRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
            b'OQGNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random
            #b'OQGRMS2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],  # buggy on isat random
            #b'OQGNO2'  : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_4],  # buggy on isat random

            # MPC Forces
            # these are NX tables
            # OQGM1 - mpc forces in the nodal frame
            b'OQMG1'   : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_forces],
            b'OQMATO1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_ato],
            b'OQMCRM1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_crm],
            b'OQMPSD1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_psd],
            b'OQMRMS1' : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_rms],
            b'OQMNO1'  : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_no],

            b'OQMG2'   : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_forces], # big random
            b'OQMATO2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_ato],
            b'OQMCRM2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_crm],
            b'OQMPSD2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_psd],
            b'OQMRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
            b'OQMNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random
            #b'OQMRMS2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_rms],  # buggy on isat random
            #b'OQMNO2'  : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_no],  # buggy on isat random

            b'ROQGM1'  : [reader_oqg._read_oqg1_3, reader_oqg._read_oqg_mpc_forces],  # relative OQM

            # acoustic pressure
            b'OPRATO1' : [reader_opr._read_opr1_3, reader_opr._read_opr_ato],
            b'OPRCRM1' : [reader_opr._read_opr1_3, reader_opr._read_opr_crm],
            b'OPRPSD1' : [reader_opr._read_opr1_3, reader_opr._read_opr_psd],
            b'OPRRMS1' : [self._nx_table_passer, self._table_passer],
            b'OPRNO1' : [reader_opr._read_opr1_3, reader_opr._read_opr_no],

            b'OPRATO2' : [reader_opr._read_opr2_3, reader_opr._read_opr_ato],
            b'OPRCRM2' : [reader_opr._read_opr2_3, reader_opr._read_opr_crm],
            b'OPRPSD2' : [reader_opr._read_opr2_3, reader_opr._read_opr_psd],
            b'OPRRMS2' : [self._nx_table_passer, self._table_passer],
            b'OPRNO2' : [reader_opr._read_opr2_3, reader_opr._read_opr_no],

            # stress
            b'OESATO1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESCRM1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESPSD1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESRMS1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESNO1'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],

            # OESXRM1C : Composite element RMS stresses in SORT1 format for random analysis that includes von Mises stress output.
            b'OESXRMS1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESXRM1C' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESXNO1' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],
            b'OESXNO1C' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],


            b'OESATO2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            b'OESCRM2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            b'OESPSD2' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],
            #b'OESRMS2' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # buggy on isat random
            #b'OESNO2'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # buggy on isat random
            b'OESRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
            b'OESNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random

            b'OESPSD2C' : [reader_oes._read_oes2_3, reader_oes._read_oes2_4],  # isat_random (nx)
            b'OSTPSD2C' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4], # isat_random (nx)
            #=======================

            # strain
            b'OSTRATO1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],
            b'OSTRCRM1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],
            b'OSTRPSD1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],
            b'OSTRRMS1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4], # isat_random
            b'OSTRNO1' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # isat_random

            b'OSTRATO2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],
            b'OSTRCRM2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],
            b'OSTRPSD2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],
            b'OSTRRMS2' : [self._table_passer, self._table_passer], # buggy on isat random
            b'OSTRNO2' : [self._table_passer, self._table_passer],  # buggy on isat random
            #b'OSTRRMS2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4], # buggy on isat random
            #b'OSTRNO2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],  # buggy on isat random

            b'OSTRMS1C' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4], # isat_random
            b'OSTNO1C' : [reader_oes._read_oes1_3, reader_oes._read_ostr1_4],  # isat_random


            # force
            b'OEFATO1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],
            b'OEFCRM1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],
            b'OEFPSD1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],
            b'OEFRMS1' : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],
            b'OEFNO1'  : [reader_oef._read_oef1_3, reader_oef._read_oef1_4],

            b'OEFATO2' : [reader_oef._read_oef2_3, reader_oef._read_oef2_4],
            b'OEFCRM2' : [reader_oef._read_oef2_3, reader_oef._read_oef2_4],
            b'OEFPSD2' : [reader_oef._read_oef2_3, reader_oef._read_oef2_4],
            #b'OEFRMS2' : [reader_oef._read_oef2_3, reader_oef._read_oef2_4], # buggy on isat random
        }

        table_mapper.update(table_mapper_geometry)
        table_mapper.update(table_mapper_random)
        table_mapper.update(nasa95_tables)
        table_mapper.update(satk_tables)
        if self.is_nx and 0:  # pragma: no cover
            _table_mapper = {
                #b'OUGRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
                #b'OUGNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random
                b'OUGRMS2' : [reader_oug._read_oug2_3, reader_oug._read_oug_rms],  # buggy on isat random
                b'OUGNO2'  : [reader_oug._read_oug2_3, reader_oug._read_oug_no],  # buggy on isat random

                #b'OQMRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
                #b'OQMNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random
                b'OQMRMS2' : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_rms],  # buggy on isat random
                b'OQMNO2'  : [reader_oqg._read_oqg2_3, reader_oqg._read_oqg_mpc_no],  # buggy on isat random

                #b'OSTRRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
                #b'OSTRNO2' : [self._table_passer, self._table_passer],  # buggy on isat random
                b'OSTRRMS2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],  # buggy on isat random
                b'OSTRNO2' : [reader_oes._read_oes2_3, reader_oes._read_ostr2_4],  # buggy on isat random

                b'OESRMS2' : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # buggy on isat random
                b'OESNO2'  : [reader_oes._read_oes1_3, reader_oes._read_oes1_4],  # buggy on isat random
                #b'OESRMS2' : [self._table_passer, self._table_passer],  # buggy on isat random
                #b'OESNO2'  : [self._table_passer, self._table_passer],  # buggy on isat random

                b'OEFNO2'  : [reader_oef._read_oef2_3, reader_oef._read_oef2_4],
                #b'OEFNO2' : [self._table_passer, self._table_passer], # buggy on isat_random_steve2.op2
            }
            for key, value in _table_mapper.items():
                table_mapper[key] = value
            #table_mapper.update(table_mapper2)
        return table_mapper

    def _read_mpf_3(self, data, ndata: int) -> int:
        """reads table 3 (the header table)

        OFMPF2E Table of fluid mode participation factors by excitation frequencies.
        OFMPF2M Table of fluid mode participation factors by normal mode.
        OSMPF2E Table of structure mode participation factors by excitation frequencies.
        OSMPF2M Table of structure mode participation factors by normal mode.
        OPMPF2E Table of panel mode participation factors by excitation frequencies.
        OPMPF2M Table of panel mode participation factors by normal mode.
        OLMPF2E Table of load mode participation factors by excitation frequencies.
        OLMPF2M Table of load mode participation factors by normal mode.
        OGMPF2E Table of grid mode participation factors by excitation frequencies.
        OGMPF2M Table of grid mode participation factors by normal mode.
        """
        #self._set_times_dtype()
        self.nonlinear_factor = np.nan
        self.is_table_1 = True
        self.is_table_2 = False
        unused_three = self.parse_approach_code(data)
        self.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label']

        ## random code
        self.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)

        ## format code
        self.format_code = self.add_data_parameter(data, 'format_code', b'i', 9, False)

        ## number of words per entry in record
        self.num_wide = self.add_data_parameter(data, 'num_wide', b'i', 10, False)

        ## acoustic pressure flag
        self.acoustic_flag = self.add_data_parameter(data, 'acoustic_flag', b'i', 13, False)

        ## thermal flag; 1 for heat transfer, 0 otherwise
        self.thermal = self.add_data_parameter(data, 'thermal', b'i', 23, False)

        #if self.analysis_code == 1:   # statics / displacement / heat flux
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5, False)
            #self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
            #self.setNullNonlinearFactor()
        #elif self.analysis_code == 2:  # real eigenvalues
            ## mode number
            #self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## eigenvalue
            #self.eign = self.add_data_parameter(data, 'eign', b'f', 6, False)
            ## mode or cycle .. todo:: confused on the type - F1???
            #self.mode_cycle = self.add_data_parameter(data, 'mode_cycle', b'i', 7, False)
            #self.update_mode_cycle('mode_cycle')
            #self.data_names = self.apply_data_code_value('data_names', ['mode', 'eign', 'mode_cycle'])
        #elif self.analysis_code == 3: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
            #self.data_code['lsdvmn'] = self.lsdvmn
        #elif self.analysis_code == 4: # differential stiffness
            #self.lsdvmn = self.get_values(data, b'i', 5) ## load set number
        if self.analysis_code == 5:   # frequency
            # frequency
            self.node_id = self.add_data_parameter(data, 'node_id', b'i', 5, fix_device_code=True)
            self.data_names = self.apply_data_code_value('data_names', ['node_id'])
            #self.freq = self.add_data_parameter(data, 'freq', b'f', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['freq'])
        #elif self.analysis_code == 6:  # transient
            ## time step
            #self.dt = self.add_data_parameter(data, 'dt', b'f', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['dt'])
        #elif self.analysis_code == 7:  # pre-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        #elif self.analysis_code == 8:  # post-buckling
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            ## real eigenvalue
            #self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            #self.data_names = self.apply_data_code_value('data_names', ['lsdvmn', 'eigr'])
        #elif self.analysis_code == 9:  # complex eigenvalues
            ## mode number
            #self.mode = self.add_data_parameter(data, 'mode', b'i', 5)
            ## real eigenvalue
            #self.eigr = self.add_data_parameter(data, 'eigr', b'f', 6, False)
            ## imaginary eigenvalue
            #self.eigi = self.add_data_parameter(data, 'eigi', b'f', 7, False)
            #self.data_names = self.apply_data_code_value('data_names', ['mode', 'eigr', 'eigi'])
        #elif self.analysis_code == 10:  # nonlinear statics
            ## load step
            #self.lftsfq = self.add_data_parameter(data, 'lftsfq', b'f', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['lftsfq'])
        #elif self.analysis_code == 11:  # old geometric nonlinear statics
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        #elif self.analysis_code == 12:  # contran ? (may appear as aCode=6)  --> straight from DMAP...grrr...
            ## load set number
            #self.lsdvmn = self.add_data_parameter(data, 'lsdvmn', b'i', 5)
            #self.data_names = self.apply_data_code_value('data_names', ['lsdvmn'])
        else:
            msg = f'invalid analysis_code...analysis_code={self.analysis_code}\ndata={self.data_code}'
            raise RuntimeError(msg)

        #print self.code_information()
        #
        self.fix_format_code()
        if self.num_wide == 8:
            self.format_code = 1
            self.data_code['format_code'] = 1
        else:
            #self.fix_format_code()
            if self.format_code == 1:
                self.format_code = 2
                self.data_code['format_code'] = 2
            assert self.format_code in [2, 3], self.code_information()

        self._parse_thermal_code()
        if self.is_debug_file:
            self.binary_debug.write('  approach_code  = %r\n' % self.approach_code)
            self.binary_debug.write('  tCode          = %r\n' % self.tCode)
            self.binary_debug.write('  isubcase       = %r\n' % self.isubcase)
        self._read_title(data)
        self._write_debug_bits()

    def _read_mpf_4(self, data: bytes, ndata: int):
        """unused"""
        if self.read_mode == 1: # or self.table_name_str not in ['OFMPF2M']:
            return ndata
        #print(self.table_name_str, ndata, self.num_wide)  # 176
        #self.show_ndata(100, types='ifs')

        structi = Struct('fiff')
        nelements = ndata // 16
        ndev = ndata % 16
        assert ndev == 0, ndev

        for i in range(nelements):
            datai = data[i*16 : (i+1)*16]
            freq, dunno_int, mag, phase = structi.unpack(datai)
            assert dunno_int == 2, str(self.node_id, freq, dunno_int, mag, phase)
            #print(self.node_id, freq, dunno_int, mag, phase)
        #print()
        if self.isubtable == -4:
            self.log.warning('%s results were read, but not saved' % self.table_name_str)
        return ndata

    def _read_pvto_3(self, data: bytes, ndata: int):
        """unused"""
        raise RuntimeError(self.read_mode)

    def _read_pvto_4(self, data: bytes, ndata: int) -> int:
        """
        Reads PARAM cards

        data = (
            AUTOSPC, 3, YES,
            GRDPNT, 1, 0,
            K6ROT, 2, 100.0,
            OUGCORD, 3, GLOBAL,
            POST, 1, -1,
            POST, 1, -2,
            POSTEXT, 3, YES,
            K6ROT, 2, 100.0,
            OIBULK, 3, YES,
            OMACHPR, 3, YES,
            POSTEXT, 3, YES,
            UNITSYS, 3, MN-MM)
        """
        if self.read_mode == 2:
            return ndata

        iloc = self.f.tell()
        try:
            ndata2 = self._read_pvto_4_helper(data, ndata)
        except (NotImplementedError, AssertionError) as error:
            #raise  # only for testing
            if 'dev' in __version__ and self.IS_TESTING or 1:
                raise  # only for testing
            self.log.error(str(error))
            log_exc(self.log)
            self.f.seek(iloc)
            ndata2 = ndata

        if 'NXVER' in self.params and not self.is_nx:
            self.set_as_nx()
            self.log.debug('found PARAM,NXVER -> setting as NX')
        return ndata2

    def _read_pvto_4_helper(self, data: bytes, ndata: int) -> int:
        """reads PARAM cards"""
        xword = (4 * self.factor)
        nvalues = ndata // xword
        assert ndata % xword == 0, ndata

        if self.size == 4:
            structi = self.struct_i
            structf = Struct(b'f')
            structs8 = self.struct_8s
            #struct2s8 = Struct(b'4s8s')
            #struct2i = self.struct_2i
            struct2f = Struct(b'ff')
            #struct2d = Struct(b'dd')
        else:
            structi = self.struct_q
            structf = Struct(b'd')
            #struct2i = self.struct_2q
            structs8 = self.struct_16s
            struct2f = Struct(b'dd')

        i = 0
        #print('---------------------------')
        #self.show_data(data, types='ifsqL')
        while i < nvalues:
            #print('-----------------------------------------------------------')
            #print('*i=%s nvalues=%s' % (i, nvalues))
            istart = i*xword
            #self.show_data(data[istart:istart+32], types='sqd')
            #self.show_data(data[istart:istart+64], types='sqd')
            if self.size == 4:
                word = data[istart:(i+2)*xword].rstrip()
            elif self.size == 8:
                bword = data[istart:(i+2)*xword]
                word = reshape_bytes_block(bword).rstrip()
            else:
                raise RuntimeError(self.size)
            key = word.decode('latin1')
            flag_data = data[(i+2)*xword:(i+3)*xword]
            flag = structi.unpack(flag_data)[0]

            # =1
            # 4 INT     I* =2
            # 4 REAL    RS* =3
            # 4 BCD(2)  CHAR4* =4
            # 4 REALDBL RD
            # =5 4 CMPLXS  CS*
            # =6 4 CMPLXD  CD
            # =7 4 LOGICAL LOGIC*

            # ----------------
            #wrong...
            # 1 PARAM(2) CHAR4
            # 3 I
            # =1 4 INT     I*
            # =2 4 REAL    RS*
            # =3 4 BCD(2)  CHAR4*
            # =4 4 REALDBL RD
            # =5 4 CMPLXS  CS*
            # =6 4 CMPLXD  CD
            # =7 4 LOGICAL LOGIC*

            #print(f'word={word!r} flag={flag}')
            #word = s8.unpack(word)[0]#.decode(self._encoding)

            #if flag == 1:
                #flag_str = 'int'
            #elif flag == 2:
                #flag_str = 'float'
            #elif flag == 3:
                #flag_str = 'str'

            # the first two entries are typically trash, then we can get values
            if flag == 1: # int
                #self.show_data(data[i*xword:(i+4)*xword], types='isq', endian=None, force=False)
                assert self.size in [4, 8], (key, self.size, flag)
                #assert word in INT_PARAMS_1, f'word={word}'
                slot = data[(i+3)*xword:(i+4)*xword]
                i += 4
                #slot = data[(i+4)*xword:(i+5)*xword]
                #i += 5
                value, = structi.unpack(slot)
                values = [value]
            elif flag == 2: # float
                assert self.size in [4, 8], (key, self.size, flag)
                slot = data[(i+3)*xword:(i+4)*xword]
                value, = structf.unpack(slot)
                values = [value]
                #assert word in FLOAT_PARAMS_1, f'word={word}'
                i += 4

            elif flag == 3: # float / string
                assert self.size in [4, 8], (key, self.size, flag)
                #slot = data[(i+3)*xword:(i+4)*xword]
                #i += 4
                slot = data[(i+3)*xword:(i+5)*xword]
                try:
                    bvalue, = structs8.unpack(slot)
                    if self.size == 8:
                        bvalue = reshape_bytes_block(bvalue)
                    value = bvalue.decode('latin1').rstrip()
                    if value:
                        if word == b'NXVER':
                            assert value.replace('.', '').isalnum(), f'{key} = {value!r}'
                        elif word == b'UNITSYS':
                            assert value.replace('-', '').isalnum(), f'{key} = {value!r}'
                        else:
                            assert value.isalnum(), f'{key} = {value!r}'
                except AssertionError:
                    value, = structf.unpack(slot[4:])
                values = [value]

                if isinstance(value, str):
                    assert word in STR_PARAMS_1, f'word={word}'
                else:
                    #if self.size == 4:
                        #self.show_data(data[istart:istart+20], types='sifqd')
                    #elif self.size == 8:
                        #self.show_data(data[istart:istart+40], types='sifqd')
                    assert word in FLOAT_PARAMS_1, f'float/str; word={word} value={value}'
                i += 5
            #elif flag == 3: # string
                #assert self.size in [4, 8], (key, self.size, flag)
                #slot = data[(i+3)*xword:(i+5)*xword]
                ##self.show_data(slot)
                #assert word in STR_PARAMS_1, f'word={word}'
                #i += 5

            elif flag == 5:  # CMPLXS  CS - FLOAT_PARAMS_2
                assert self.size in [4, 8], (key, self.size, flag)
                slot = data[(i+3)*xword:(i+5)*xword]
                #self.show_data(data[(i+3)*xword:(i+5)*xword], types='ifsqd', endian=None, force=False)
                values = struct2f.unpack(slot)
                values = list(values)
                assert word in FLOAT_PARAMS_2, f'word={word}'
                i += 5

            elif flag == 7: # logical/int
                assert self.size in [4, 8], (key, self.size, flag)
                slot = data[(i+3)*xword:(i+4)*xword]
                value, = structi.unpack(slot)
                values = [value]
                i += 4
            else:
                self.show_data(data[i*xword:], types='ifsqd', endian=None, force=False)
                self.log.error('%r' % word)
                raise NotImplementedError(f'{word!r} is not a supported PARAM; flag={flag}')

            #i, value = self._old_pvto(word, data, i, xword,
                                      #struct2i, struct2f, structs8)

            param = PARAM(key, values, comment='')
            self.params[key] = param
            #print(f'{key} ({flag}) = {value!r}')
            del key, values
            #print(param.rstrip())
        return nvalues

    def _not_available(self, data: bytes, ndata: int):
        """testing function"""
        if ndata > 0:
            raise RuntimeError('this should never be called...'
                               'table_name={self.table_name!r} len(data)={ndata}')

    def _table_crasher(self, data: bytes, ndata: int):
        """auto-table crasher"""
        if self.is_debug_file:
            self.binary_debug.write(f'  crashing table = {self.table_name}\n')
            raise NotImplementedError(self.table_name)
        return ndata

    def _nx_table_passer(self, data, ndata: int):
        """auto-table skipper"""
        self.to_nx(f' because table_name={self.table_name} was found')
        self._table_passer(data, ndata)

    def _table_passer(self, data, ndata: int):
        """auto-table skipper"""
        if self.is_debug_file:
            self.binary_debug.write(f'  skipping table = {self.table_name}\n')
        if self.table_name not in GEOM_TABLES and self.isubtable > -4:
            self.log.warning(f'    skipping table: {self.table_name_str}')
        if not is_release and self.isubtable > -4:
            if self.table_name in GEOM_TABLES and not self.make_geom:
                pass
            else:
                print(f'dont skip table {self.table_name_str!r}')
                raise RuntimeError(f'dont skip table {self.table_name_str!r}')
        return ndata

    def _validate_op2_filename(self, op2_filename: Optional[str]) -> str:
        """
        Pops a GUI if the op2_filename hasn't been set.

        Parameters
        ----------
        op2_filename : str
            the filename to check (None -> gui)

        Returns
        -------
        op2_filename : str
            a valid file string

        """
        if op2_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran OP2 (*.op2)|*.op2|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran OP2 (*.op2);;All files (*)"
            title = 'Please select a OP2 to load'
            op2_filename, unused_wildcard_level = load_file_dialog(
                title, wildcard_wx, wildcard_qt, dirname='')
            assert op2_filename is not None, op2_filename
        return op2_filename

    def _create_binary_debug(self):
        """Instatiates the ``self.binary_debug`` variable/file"""
        if hasattr(self, 'binary_debug') and self.binary_debug is not None:
            self.binary_debug.close()
            del self.binary_debug

        self.is_debug_file, self.binary_debug = create_binary_debug(
            self.op2_filename, self.debug_file, self.log)

    def _setup_filenames(self, op2_filename: Optional[str], force: bool=True):
        if op2_filename:
            fname = os.path.splitext(op2_filename)[0]
            self.op2_filename = op2_filename
            self.bdf_filename = fname + '.bdf'
            self.f06_filename = fname + '.f06'
            self.des_filename = fname + '.des'
            self.h5_filename = fname + '.h5'

    def read_op2(self, op2_filename=None,
                 combine: bool=False,
                 load_as_h5: bool=False,
                 h5_file=None,
                 mode: Optional[str]=None) -> None:
        """
        Starts the OP2 file reading

        Parameters
        ----------
        op2_filename : str
            the op2 file
        combine : bool; default=True
            True : objects are isubcase based
            False : objects are (isubcase, subtitle) based;
                    will be used for superelements regardless of the option
        load_as_h5 : default=None
            False : don't setup the h5_file
            True : loads the op2 as an h5 file to save memory
                   stores the result.element/data attributes in h5 format
        h5_file : h5File; default=None
            None : ???
            h5File : ???

        +--------------+-----------------------+
        | op2_filename | Description           |
        +--------------+-----------------------+
        |     None     | a dialog is popped up |
        +--------------+-----------------------+
        |    string    | the path is used      |
        +--------------+-----------------------+
        """
        self._setup_filenames(op2_filename, force=False)
        self.op2_reader.load_as_h5 = load_as_h5
        if load_as_h5:
            h5_file = None
            import h5py
            self.h5_file = h5py.File(self.h5_filename, 'w')
            self.op2_reader.h5_file = self.h5_file

        self._count = 0
        if self.read_mode == 1:
            #sr = list(self._results.saved)
            #sr.sort()
            #self.log.debug('_results.saved = %s' % str(sr))
            #self.log.info('_results.saved = %s' % str(sr))
            pass

        if self.read_mode != 2:
            op2_filename = self._validate_op2_filename(op2_filename)
            self.log.info(f'op2_filename = {op2_filename!r}')
            self._setup_filenames(op2_filename, force=True)
            if not is_binary_file(op2_filename):
                if os.path.getsize(op2_filename) == 0:
                    raise IOError(f'op2_filename={op2_filename!r} is empty.')
                raise IOError(f'op2_filename={op2_filename!r} is not a binary OP2.')

        self._create_binary_debug()
        self._setup_op2()
        _op2 = self.op2_reader.op2
        #is_nasa_nastran = False
        #if is_nasa_nastran:
            #self.show(104, types='ifs', endian=None)
            #self.show(52, types='ifs', endian=None)
            #aa
            #data = _op2.f.read(4)
            #_op2.n += 8
            #_op2.f.seek(_op2.n)
        #else:
        self.op2_reader.read_nastran_version(mode)
        data = _op2.f.read(4)
        _op2.f.seek(_op2.n)
        if len(data) == 0:
            raise FatalError('There was a Nastran FATAL Error.  Check the F06.\n'
                             'No tables exist...check for a license issue')

        #=================
        table_name = self.op2_reader._read_table_name(rewind=True, stop_on_failure=False)
        if table_name is None:
            raise FatalError('There was a Nastran FATAL Error.  Check the F06.\n'
                             'No tables exist...check for a license issue')

        self._make_tables()
        table_names = []
        try:
            self._read_tables(table_name, table_names)
        except EmptyRecordError:
            self.show(500, types='ifs', endian=None, force=False)
            raise
            op2_reader = self.op2_reader
            op2_reader.read_markers([1, 0, 0, 0])
            self.show(500, types='ifs', endian=None, force=False)
            self._finish()

        self.close_op2(force=False)
        #self.remove_unpickable_data()
        return table_names

    def close_op2(self, force=True):
        """closes the OP2 and debug file"""
        if self.is_debug_file:
            self.binary_debug.write('-' * 80 + '\n')
            self.binary_debug.write('f.tell()=%s\ndone...\n' % self.f.tell())
            self.binary_debug.close()

        if self._close_op2 or force:
            if self.f is not None:
                # can happen if:
                #  - is ascii file
                self.f.close()
            del self.binary_debug
            del self.f
            self._cleanup_data_members()
            self._cleanup_words()
            #self.op2_reader.h5_file.close()

    def _cleanup_words(self):
        """
        Remove internal parameters that are not useful and just clutter
        the object attributes.
        """
        words = [
            'isubcase', 'int3', '_table4_count', 'nonlinear_factor',
            'is_start_of_subtable', 'superelement_adaptivity_index',
            'thermal_bits', 'is_vectorized', 'pval_step', #'_frequencies',
            '_analysis_code_fmt', 'isubtable', '_data_factor', 'sort_method',
            'acoustic_flag', 'approach_code', 'format_code_original',
            'element_name', 'sort_bits', 'code', 'n', 'use_vector', 'ask',
            'stress_bits', 'expected_times', 'table_code', 'sort_code',
            'is_all_subcases', 'num_wide', '_table_mapper', 'label',
            'apply_symmetry',
            'words', 'device_code', 'table_name', '_count', 'additional_matrices',
            # 350
            'data_names', '_close_op2',
            'op2_reader',
            # 74
            'generalized_tables',
            # 124
            'is_table_1', 'is_table_2', 'ntotal', 'element_mapper',
            'is_debug_file', 'debug_file',
            '_results', 'skip_undefined_matrices',
            # 140
            #---------------------------------------------------------
            # dont remove...
            # make_geom, title, read_mode
            # result_names, op2_results

        ]
        for word in words:
            if hasattr(self, word):
                delattr(self, word)

    def _setup_op2(self):
        """
        Does preliminary op2 tasks like:
          - open the file
          - set the endian
          - preallocate some struct objects

        """
        #: file index
        self.n = 0
        self.table_name = None


        if not hasattr(self, 'f') or self.f is None:
            #: the OP2 file object
            op2_filename = self.op2_filename
            self.f = open(op2_filename, 'rb')
            #: the endian in bytes
            self._endian = None
            #: the endian in unicode
            self._uendian = None
            flag_data = self.f.read(20)
            self.f.seek(0)

            #(4, 3, 4, 24) ???
            #(8, 3, 0, 8, 24)
            little_data = unpack(b'<5i', flag_data)
            big_data = unpack(b'>5i', flag_data)
            if big_data[0] in [4, 8]:
                self._uendian = '>'
                self._endian = b'>'
                size = big_data[0]
            elif little_data[0] in [4, 8]:
                self._uendian = '<'
                self._endian = b'<'
                size = little_data[0]
            #elif unpack(b'<ii', flag_data)[0] == 4:
                #self._endian = b'<'
            else:
                # Matrices from test show
                # (24, 10, 10, 6, 2) before the Matrix Name...
                print(little_data, big_data)
                self.show(30, types='ifs', endian='<')
                self.show(30, types='ifs', endian='>')
                self.show(12, types='ifs', endian='<')
                self.show(12, types='ifs', endian='>')
                #self.show_data(flag_data, types='iqlfsld', endian='<')
                #print('----------')
                #self.show_data(flag_data, types='iqlfsld', endian='>')
                raise FatalError('cannot determine endian')
        else:
            self.op2_reader._goto(self.n)

        if self.read_mode == 1:
            self._set_structs(size)

    def _make_tables(self):
        return
        #global RESULT_TABLES, NX_RESULT_TABLES, MSC_RESULT_TABLES
        #table_mapper = self._get_table_mapper()
        #RESULT_TABLES = table_mapper.keys()

    def _read_tables(self, table_name: bytes, table_names: List[bytes]) -> None:
        """
        Reads all the geometry/result tables.
        The OP2 header is not read by this function.

        Parameters
        ----------
        table_name : bytes str
            the first table's name
        table_names : List[bytes str]
            the table names that were read

        """
        op2_reader = self.op2_reader
        self.table_count = defaultdict(int)
        while table_name is not None:
            self.table_count[table_name] += 1
            table_names.append(table_name)

            if self.is_debug_file:
                self.binary_debug.write('-' * 80 + '\n')
                self.binary_debug.write(f'table_name = {table_name!r}\n')

            if is_release:
                self.log.debug(f'  table_name={table_name!r}')

            self.table_name = table_name
            #if 0:
                #op2_reader._skip_table(table_name)
            #else:
            #print(table_name, table_name in op2_reader.mapped_tables)
            if table_name in self.generalized_tables:
                t0 = self.f.tell()
                self.generalized_tables[table_name](self)
                assert self.f.tell() != t0, 'the position was unchanged...'
            elif table_name in op2_reader.mapped_tables:
                t0 = self.f.tell()
                op2_reader.mapped_tables[table_name]()
                assert self.f.tell() != t0, 'the position was unchanged...'
            elif table_name in GEOM_TABLES:
                op2_reader.read_geom_table()  # DIT (agard)
            elif table_name in MATRIX_TABLES:
                op2_reader.read_matrix(table_name)
            elif table_name in RESULT_TABLES:
                op2_reader.read_results_table()
            elif self.skip_undefined_matrices:
                op2_reader.read_matrix(table_name)
            elif table_name.strip() in self.additional_matrices:
                op2_reader.read_matrix(table_name)
            else:
                #self.show(1000, types='ifsq')
                msg = (
                    f'Invalid Table = {table_name!r}\n\n'
                    'If you have matrices that you want to read, see:\n'
                    '  model.set_additional_matrices_to_read(matrices)\n'
                    '  matrices = {\n'
                    "      b'BHH' : True,\n"
                    "      b'KHH' : False,\n"
                    '  }  # you want to read some matrices, but not others\n'
                    "  matrices = [b'BHH', b'KHH']  # assumes True\n\n"

                    'If you the table is a geom/result table, see:\n'
                    '  model.set_additional_result_tables_to_read(methods_dict)\n'
                    "  methods_dict = {\n"
                    "      b'OUGV1' : [method3, method4],\n"
                    "      b'GEOM4SX' : [method3, method4],\n"
                    "      b'OES1X1' : False,\n"
                    '  }\n\n'

                    'If you want to take control of the OP2 reader (mainly useful '
                    'for obscure tables), see:\n'
                    "  methods_dict = {\n"
                    "      b'OUGV1' : [method],\n"
                    '  }\n'
                    '  model.set_additional_generalized_tables_to_read(methods_dict)\n'
                )
                raise NotImplementedError(msg)

            table_name = op2_reader._read_table_name(last_table_name=table_name,
                                                     rewind=True, stop_on_failure=False)

    def set_additional_generalized_tables_to_read(self, tables):
        """
        Adds methods to call a generalized table.
        Everything is left to the user.

        ::

          def read_some_table(self):
              # read the data from self.f
              pass

          # let's overwrite the existing OP2 table
          model2 = OP2Geom(debug=True)
          generalized_tables = {
              b'GEOM1S' : read_some_table,
          }

          model.set_additional_generalized_tables_to_read(generalized_tables)

        """
        self._update_generalized_tables(tables)
        self.generalized_tables = tables

    def set_additional_result_tables_to_read(self, tables):
        """
        Adds methods to read additional result tables.
        This is expected to really only be used for skipping
        unsupported tables or disabling enabled tables that are
        buggy (e.g., OUGV1).

        Parameters
        ----------
        tables : Dict[bytes] = varies
            a dictionary of key=name, value=list[method3, method4]/False,
            False : skips a table
                applies self._table_passer to method3 and method4
            method3 : function
                function to read table 3 results (e.g., metadata)
            method4 : function
                function to read table 4 results (e.g., the actual results)

        """
        self._update_generalized_tables(tables)
        table_mapper = self._get_table_mapper()
        #is_added = False
        def func():
            """overloaded version of _get_table_mapper"""
            #if is_added:
                #return table_mapper
            for _key, methods in tables.items():
                if methods is False:
                    table_mapper[_key] = [self._table_passer, self._table_passer]
                else:
                    assert len(methods) == 2, methods
                    table_mapper[_key] = methods
            #is_added = True
            return table_mapper
        self._get_table_mapper = func

    def _update_generalized_tables(self, tables):
        """
        helper function for:
         - set_additional_generalized_tables_to_read
         - set_additional_result_tables_to_read

        """
        global NX_RESULT_TABLES
        global MSC_RESULT_TABLES
        global RESULT_TABLES
        failed_keys = []
        keys = list(tables.keys())
        for _key in keys:
            if not isinstance(_key, bytes):
                failed_keys.append(_key)
            if hasattr(self, 'is_nx') and self.is_nx:
                NX_RESULT_TABLES.append(_key)
            else:
                MSC_RESULT_TABLES.append(_key)
        if failed_keys:
            failed_keys_str = [str(_key) for _key in failed_keys]
            raise TypeError('[%s] must be bytes' % ', '. join(failed_keys_str))
        RESULT_TABLES = NX_RESULT_TABLES + MSC_RESULT_TABLES

        #RESULT_TABLES.sort()
        #assert 'OESXRMS1' in RESULT_TABLES, RESULT_TABLES

    def set_additional_matrices_to_read(self, matrices: Union[List[str], Dict[str, bool]]):
        """
        Matrices (e.g., KHH) can be sparse or dense.

        Parameters
        ----------
        matrices : List[str]; Dict[str] = bool
            List[str]:
                simplified method to add matrices; value will be True
            Dict[str] = bool:
                a dictionary of key=name, value=True/False,
                where True/False indicates the matrix should be read

        .. note:: If you use an already defined table (e.g. KHH), it
                  will be ignored.  If the table you requested doesn't
                  exist, there will be no effect.
        .. note:: Do not use this for result tables like OUGV1, which
                  store results like displacement.  Those are not matrices.
                  Matrices are things like DMIGs.

        """
        if isinstance(matrices, list):
            matrices2 = {}
            for matrix in matrices:
                assert isinstance(matrix, str), 'matrix=%r' % str(matrix)
                matrices2[matrix] = True
            matrices = matrices2

        self.additional_matrices = matrices
        self.additional_matrices = {}
        for matrix_name, matrix in matrices.items():
            if isinstance(matrix_name, bytes):
                self.additional_matrices[matrix_name] = matrix
            else:
                self.additional_matrices[matrix_name.encode('latin1')] = matrix

    def _finish(self):
        """
        Clears out the data members contained within the self.words variable.
        This prevents mixups when working on the next table, but otherwise
        has no effect.

        """
        for word in self.words:
            if word != '???' and hasattr(self, word):
                if word not in ['Title', 'reference_point']:
                    delattr(self, word)
        self.obj = None
        if hasattr(self, 'subtable_name'):
            del self.subtable_name

    def _read_psdf_3(self, data: bytes, ndata: int):
        """reads the PSDF table"""
        #(50, 2011, 4001, 0, 302130, 3
        # strip off the title
        unused_three = self.parse_approach_code(data)
        self.words = [
            'approach_code', 'table_code', '???', 'isubcase',
            '???', '???', '???', 'random_code',
            'format_code', 'num_wide', '???', '???',
            'acoustic_flag', '???', '???', '???',
            '???', '???', '???', '???',
            '???', '???', 'thermal', '???',
            '???', 'Title', 'subtitle', 'label'
        ]

        ## random code
        self.random_code = self.add_data_parameter(data, 'random_code', b'i', 8, False)
        self._read_title(data)

        # simplifying to see the data better
        del self.data_code['title']
        del self.data_code['label']
        del self.data_code['subtitle']
        del self.data_code['subtitle_original']
        del self.data_code['superelement_adaptivity_index']
        #del self.data_code['pval_step']
        del self.data_code['table_name']

        del self.data_code['_encoding']
        del self.data_code['load_as_h5']
        del self.data_code['h5_file']
        del self.data_code['is_msc']
        #del self.data_code['is_nasa95']
        del self.data_code['pval_step']

        # wrong
        del self.data_code['isubcase']
        #del self.data_code['random_code']
        #del self.data_code['sort_bits']
        #del self.data_code['device_code']
        #del self.data_code['sort_code']
        #del self.data_code['sort_method']
        #print(self.data_code)

        if self.read_mode == 1:
            return ndata
        # just stripping off title
        #self.show_data(data[:200], types='if')

        #  stripping off zeros
        #self.show_data(data[:52], types='ifs')

        #self.show_data(data[:40], types='if')

        approach_code, tcode, int3, frame_id, int5, dof, float7, rms_value, float9, int10, stress_strain_flag = unpack(
            self._endian + b'6i 3f 2i', data[:44])
        self.stress_strain_flag = stress_strain_flag

        ints = np.frombuffer(data[:200], dtype=self.idtype)
        if ints[11:].max() > 0:
            self.log.warning(f'ints11 = {ints[11:].tolist()}')

        node = int5 // 10
        #dof = int5 % 10
        #from pyNastran.op2.op2_interface.op2_codes import TABLE_CODE_MAP
        #title = self.title
        #subtitle = self.subtitle
        #label = self.label
        #approach_code={iapproach_code}  tcode={tcode} table_code={self.table_code}
        #print(f'analysis_code={self.analysis_code} '
        #print(f'title={title!r} subtitle={subtitle!r} label={label!r}')

        if (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 1, 0):
            word = 'displacements'
        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 2, 0):
            word = 'load_vectors'
        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 3, 0):
            word = 'spc_forces'
        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 4, 0):
            word = 'force'

        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 5, 0):
            word = 'stress'
        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 5, 2):
            word = 'strain'

        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 10, 0):
            word = 'velocities'
        elif (self.analysis_code, self.table_code, self.stress_strain_flag) == (5, 11, 0):
            word = 'accelerations'
        else:  # pragma: no cover
            #print(f'table_code={self.table_code} table={TABLE_CODE_MAP[self.table_code]!r}')
            print(f'analysis_code={self.analysis_code} approach_code={approach_code} tcode={tcode} table_code={self.table_code} '
                  f'int3={int3} frame_id={frame_id} node={node} dof={dof} '
                  f'float7={float7} rms_value={rms_value:.5e} float9={float9:.4e} int10={int10} stress_strain_flag={stress_strain_flag}')
            raise NotImplementedError(f'analysis_code={self.analysis_code} '
                                      f'table_code={self.table_code} '
                                      f'stress_strain_flag={self.stress_strain_flag} is not supported')

        self.node = node
        self.dof = dof
        self.word = word
        return ndata
        #self.show_data(data, types='ifs', endian=None)
        #aaaa

    def _read_psdf_4(self, data: bytes, ndata: int):
        """reads the PSDF table"""
        if self.read_mode == 1:
            return ndata
        #self.show_data(data[:100], types='ifs', endian=None)
        data2 = np.frombuffer(data, dtype=self.fdtype)
        ndata = len(data2)
        nfreqs = ndata // 2
        data2 = data2.reshape(nfreqs, 2)
        #last2 = data2[-2:, 1]
        #self.log.warning(f'skipping PSDF; nfreqs={nfreqs} [{last2[0]:.6e},{last2[1]:.6e}] '
                         #f'ymin={data2[:,1].min():.6e} ymax={data2[:,1].max():.6e}') #  {self.data_code}
        # self.show_data(), self._read_psdf_4
        key = (self.label, self.node, self.dof)
        slot = getattr(self.op2_results.psds, self.word)
        assert key not in slot, slot
        slot[key] = data2
        del self.node
        del self.dof
        del self.word

def main():  # pragma: no cover
    """testing pickling"""
    from pickle import dump, load
    txt_filename = 'solid_shell_bar.txt'
    pickle_file = open(txt_filename, 'wb')
    op2_filename = 'solid_shell_bar.op2'
    op2 = OP2_Scalar()
    op2.read_op2(op2_filename)
    #print(op2.displacements[1])
    dump(op2, pickle_file)
    pickle_file.close()

    pickle_file = open(txt_filename, 'r')
    op2 = load(pickle_file)
    pickle_file.close()
    #print(op2.displacements[1])


    #import sys
    #op2_filename = sys.argv[1]

    #o = OP2_Scalar()
    #o.read_op2(op2_filename)
    #(model, ext) = os.path.splitext(op2_filename)
    #f06_outname = model + '.test_op2.f06'
    #o.write_f06(f06_outname)

def create_binary_debug(op2_filename: str, debug_file: str, log) -> Tuple[bool, Any]:
    """helper method"""
    binary_debug = None

    if debug_file is not None:
        #: an ASCII version of the op2 (creates lots of output)
        log.debug('debug_file = %s' % debug_file)
        binary_debug = open(debug_file, 'w')
        binary_debug.write(op2_filename + '\n')
        is_debug_file = True
    else:
        is_debug_file = False
    return is_debug_file, binary_debug


if __name__ == '__main__':  # pragma: no cover
    main()
