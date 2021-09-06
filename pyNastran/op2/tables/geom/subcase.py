import numpy as np
from pyNastran.bdf.subcase import Subcase
from pyNastran.op2.op2_interface.op2_reader import reshape_bytes_block # , reshape_bytes_block_size


def set_casecc(self, data: bytes, idtype: str, fdtype: str, size: int=4,
               nastran_format: str='nx'):
    """
    Word Name Type Description
    subcase, mpc, spc, load, method_structure, deform, temp_load, temp_mat_init, tic, nlload_set, nlload_media, nlload_format, dload, freq, tfl = ints[:15]
    1 SID                I Subcase identification number
    2 MPCSET             I Multipoint constraint set (MPC)
    3 SPCSET             I Single point constraint set (SPC)
    4 ESLSET             I External static load set (LOAD)
    5 REESET             I Real eigenvalue extraction set (METHOD(STRUCTURE))
    6 ELDSET             I Element deformation set (DEFORM)
    7 THLDSET            I Thermal load set (TEMP(LOAD))
    8 THMATSET           I Thermal material set TEMP(MAT or INIT)
    9 TIC                I Transient initial conditions (IC)
    10 NONPTSET          I Nonlinear load output set (NLLOAD)
    11 NONMEDIA          I Nonlinear load output media (NLLOAD)
    12 NONFMT            I Nonlinear load output format (NLLOAD)
    13 DYMLDSET          I Dynamic load set (DLOAD)
    14 FEQRESET          I Frequency response set (FREQUENCY)
    15 TFSET             I Transfer function set (TFL)
    symflag, oload_set, oload_media, oload_format, disp_set, disp_media, disp_format, stress_set, stress_media, stress_set, force_set, force_media, force_set = ints[15:30]
    16 SYMFLG            I Symmetry flag (SYMSEQ and SUBSEQ)
    17 LDSPTSET          I Load output set (OLOAD)
    18 LDSMEDIA          I Load output media (OLOAD)
    19 LDSFMT            I Load output format (OLOAD)
    20 DPLPTSET          I Displ., temp., or pressure output set (DISP,THERM,PRES)
    21 DPLMEDIA          I Displ., temp., or pressure output media (DISP,THERM,PRES)
    22 DPLFMT            I Displ., temp., or pressure output format (DISP,THERM,PRES)
    23 STSPTSET          I Stress output set (STRESS)
    24 STSMEDIA          I Stress output media (STRESS)
    25 STSFMT            I Stress output format (STRESS)
    26 FCEPTSET          I Force (or flux) output set (FORCE or FLUX)
    27 FCEMEDIA          I Force (or flux) output media (FORCE or FLUX)
    28 FCEFMT            I Force (or flux) output format (FORCE or FLUX)

    29 ACCPTSET          I Acceleration (or enthalpy delta) output set (ACCEL or HDOT)
    30 ACCMEDIA          I Acceleration (or enthalpy delta) output media (ACCE, HDOT)
    31 ACCFMT            I Acceleration (or enthalpy delta) output format (ACCE, HDOT)
    32 VELPTSET          I Velocity (or enthalpy) output set (VELOCITY or ENTHALPY)
    33 VELMEDIA          I Velocity (or enthalpy) output media (VELOCITY) or ENTHALPY)
    34 VELFMT            I Velocity (or enthalpy) output format (VELOCITY) or ENTHALPY)
    35 FOCPTSET          I Forces of single-point constraint output set (SPCFORCE)
    36 FOCMEDIA          I Forces of single-point constraint output media (SPCFORCE)
    37 FOCFMT            I Forces of single-point constraint output format (SPCFORCE)
    38 TSTEPTRN          I Time step set for transient analysis (TSTEP)
    39 TITLE(32)     CHAR4 Title character string (TITLE)
    71 SUBTITLE(32)  CHAR4 Subtitle character string (SUBTITLE)
    103 LABEL(32)    CHAR4 LABEL character string (LABEL)
    135 STPLTFLG         I Model plot flag: set to 1 if OUTPUT(PLOT) is specified
    136 AXSYMSET         I Axisymmetric set (AXISYMMETRIC)
    137 NOHARMON         I Number of harmonics to output (HARMONICS)
    138 TSTRV            I Need definition
    139 K2PP(2)      CHAR4 Name of direct input (p-set) stiffness matrix (K2PP)
    141 M2PP(2)      CHAR4 Name of direct input (p-set) mass matrix (M2PP)
    143 B2PP(2)      CHAR4 Name of direct input (p-set) damping matrix (B2PP)
    145 OUTRESPV         I Output frequencies or times (OFREQ or OTIME)
    146 SEDR             I Data recovery superelement list (SEDR)
    147 FLDBNDY          I Fluid boundary element selection (MFLUID)
    148 CEESET           I Complex eigenvalue extraction set (CMETHOD)
    149 DAMPTBL          I Structural damping table set (SDAMP(STRUCT)
    151 SSDSET           I Solution set displacements output set (SDISP)
    152 SSDMEDIA         I Solution set displacements output media (SDISP)
    153 SSDFMT           I Solution set displacements output format (SDISP)
    154 SSVSET           I Solution set velocities output set (SVELO)
    155 SSVMEDIA         I Solution set velocities output media (SVELO)
    156 SSVFMT           I Solution set velocities output format (SVELO)
    157 SSASET           I Solution set accelerations output set (SACCE)
    158 SSAMEDIA         I Solution set accelerations output media (SACCE)
    159 SSAFMT           I Solution set accelerations output format (SACCE)
    160 NONLINLD         I Nonlinear load set in transient problems (NONLINEAR)
    161 PARTIT           I Partitioning set (PARTN)
    162 CYCLIC           I Symmetry option in cyclic symmetry (DSYM)
    163 RANDOM           I Random analysis set (RANDOM)
    164 NONPARAM         I Nonlinear static analysis control parameters (NLPARM)
    165 FLUTTER          I Flutter set (FMETHOD)
    166 LCC              I Number of words in this record up to LSEM
    167 GPFSET           I Grid point force output set (GPFORCE)
    168 GPFMEDIA         I Grid point force output media (GPFORCE)
    169 GPFFMT           I Grid point force output format (GPFORCE)
    170 ESESET           I Strain energy output set (ESE)
    171 ESEMEDIA         I Strain energy output media (ESE)
    172 ESEFMT           I Strain energy output format (ESE)
    173 ARFPTSET         I Aerodynamic force output set (AEROF)
    174 ARFMEDIA         I Aerodynamic force output media (AEROF)
    175 ARFFMT           I Aerodynamic force output format (AEROF)
    176 SEID             I Superelement ID (SUPER)
    177 LCN              I Load column number (SUPER)
    178 GUST             I Gust load selection (GUST)
    179 SEFINAL          I Final Superelement ID (SEFINAL)
    180 SEMG             I Generate matrices (K,M,B,K4) for superelement set or ID (SEMG)
    181 SEKR             I Reduce stiffness matrix (K) for superelement set or ID (SEKR)
    182 SELG             I Generate static loads for superelement set or ID (SELG)
    183 SELR             I Reduce static loads for superelement set or ID (SELR)
    184 SEEX             I Superelement set or ID to be excluded (SEEXCLUDE)
    185 K2GG(2)      CHAR4 Name of direct input (g-set) stiffness matrix (K2GG)
    187 M2GG(2)      CHAR4 Name of direct input (g-set) stiffness matrix (M2GG)
    189 B2GG(2)      CHAR4 Name of direct input (g-set) stiffness matrix (B2GG)
    191 SVSET            I Solution eigenvector output set (SVECTOR)
    192 SVMEDIA          I Solution eigenvector output media (SVECTOR)
    193 SVFMT            I Solution eigenvectors output format (SVECTOR)
    194 FLUPTSET         I Fluid pressure output set (MPRES)
    195 FLUMEDIA         I Fluid pressure output media (MPRES)
    196 FLUFMT           I Fluid pressure output format (MPRES)
    197 HOUT(3)          I Cyclic symmetry harmonic output (HOUTPUT)
    200 NOUT(3)          I Cyclic symmetry physical output (NOUTPUT)
    203 P2G(2)       CHAR4 Name of direct input (g-set) static loads matrix (P2G)
    205 LOADSET          I Sequence of static loads sets (LOADSET)
    206 SEMR             I Generate matrices (M,B,K4) for superelement set or ID (SEMG)
    207 VONMISES         I von Mises fiber (STRESS)
    208 SECMDFLG         I Superelement command existence flag
    209 GPSPTSET         I Grid point stress output set (GPSTRESS)
    210 GPSMEDIA         I Grid point stress output media (GPSTRESS)
    211 GPSFMT           I Grid point stress output format (GPSTRESS)
    212 STFSET           I Grid point stress field output set (STRFIELD)
    213 STFMEDIA         I Grid point stress field output media (STRFIELD
    214 STFFMT           I Grid point stress field output format (STRFIELD)
    215 CLOAD            I Superelement static load combination set (CLOAD)
    216 SET2ID           I Legacy design sensitivity constraint and variable set (SET2)
    217 DSAPRT           I Legacy design sensitivity analysis print option (SENSITY)
    218 DSASTORE         I Legacy design sensitivity analysis store option (SENSITY)
    219 DSAOUTPT         I Legacy design sensitivity analysis OUTPUT4 option (SENSITY)
    220 STNSET           I Strain output set (STRAIN)
    221 STNMEDIA         I Strain output media (STRAIN)
    222 STNFMT           I Strain output format (STRAIN)
    223 APRESS           I Aerodynamic pressure output set (APRESSURE)
    224 TRIM             I Aerostatic trim variable constrain set (TRIM)
    225 MODLIST          I Output modes list set (OMODES)
    226 REESETF          I Real eigenvalue extraction set for fluid (METHOD(FLUID))
    227 ESDPTSET         I Element stress discontinuity output set (ELSDCON)
    228 ESDMEDIA         I Element stress discontinuity output media (ELSDCON)
    229 ESDFMT           I Element stress discontinuity output format (ELSDCON)
    230 GSDPTSET         I Grid point stress discontinuity output set (GPSDCON)
    231 GSDMEDIA         I Grid point stress discontinuity output media (GPSDCON)
    232 GSDFMT           I Grid point stress discontinuity output format (GPSDCON)
    233 SEDV             I Generate pseudo-loads for superelement set or identification number (SEDV)
    234 SERE             I Generate responses for superelement set or ID (SERESP)
    235 SERS             I Restart processing for superelement set or ID (SERS)
    236 CNTSET           I Slideline contact output set (BOUTPUT)
    237 CNTMEDIA         I Slideline contact output media (BOUTPUT)
    238 CNTFMT           I Slideline contact output format (BOUTPUT)
    239 DIVERG           I Aerostatic divergence control parameter set (DIVERG)
    240 OUTRCV           I P-element output control parameters (OUTRCV)
    241 STATSUBP         I Static subcase identification number for pre-load (STATSUB(PRELOAD))
    242 MODESELS         I Mode selection set identification number for the structure (MODESELECT)
    243 MODESELF         I Mode selection set identification number for the fluid (MODESELECT)
    244 SOLNUM           I Solution sequence number
    245 ANLOPT           I SOL 601 analysis option: 106, 129, 153 or 159
    246 ADAPT            I P-element adaptivity control parameter set (ADAPT)
    247 DESOBJ           I Design objective set (DESOBJ)
    248 DESSUB           I Design constraint set for current subcase (DESSUB)
    249 SUBSPAN          I DRSPAN defined set ID of DRESP1 responses specific to this subcase.
    250 DESGLB           I Design constraint set for all subcases (DESGLB)
    251 ANALYSIS     CHAR4 Type of analysis (ANALYSIS)
    252 GPQSTRS          I CQUAD4 grid point corner stress option (STRESS)
    253 GPQFORC          I CQUAD4 grid point corner force option (STRESS)
    254 GPQSTRN          I CQUAD4 grid point corner strain option (STRESS)
    255 SUPORT1          I Supported degree-of-freedom set (SUPORT1)
    256 STATSUBB         I Static subcase ID for buckling (STATSUB(BUCKLE))
    257 BCID             I Boundary condition ID (BC)
    258 AUXMODEL         I Auxiliary model ID (AUXMODEL)
    259 ADACT            I P-element adaptivity active subcase flag (ADACT)
    260 DATSET           I P-element output set (DATAREC)
    261 DATMEDIA         I P-element output media (DATAREC)
    262 DATFMT           I P-element output format (DATAREC)
    263 VUGSET           I View-grid and element output set (VUGRID)
    264 VUGMEDIA         I View-grid and element output media (VUGRID)
    265 VUGFMT           I View-grid and element output format (VUGRID)
    266 MPCFSET          I Forces of multipoint constraint output set (MPCFORCE)
    267 MPCMEDIA         I Forces of multipoint constraint output media (MPCFORCE)
    268 MPCFFMT          I Forces of multipoint constraint output format (MPCFORCE)
    269 REUESET          I Real unsymmetric eigenvalue extraction set (UMETHOD)
    270 DAMPTBLF         I Structural damping table set for the fluid (SDAMP(FLUID)
    271 ITERMETH         I Iterative solver control parameters (SMETHOD)
    272 NLSSET           I Nonlinear stress output set (NLSTRESS)
    273 NLSMEDIA         I Nonlinear stress output media (NLSTRESS)
    274 NLSFMT           I Nonlinear stress output format (NLSTRESS)
    275 MODTRKID         I Mode tracking control parameter set (MODTRAK)
    276 DSAFORM          I Design sensitivity output format: 1=yes,2=no (DSAPRT)
    277 DSAEXPO          I Design sensitivity output export: 1=no,2=yes (DSAPRT)
    278 DSABEGIN         I Design sensitivity output start iteration (DSAPRT)
    279 DSAINTVL         I Design sensitivity output interval (DSAPRT)
    280 DSAFINAL         I Design sensitivity output final iteration (DSAPRT)
    281 DSASETID         I Design sensitivity output set (DSAPRT)
    282 SORTFLG          I Overall SORT1/SORT2 flag: 1 means SORT1 and 2 means SORT2.
    283 RANDBIT          I Random analysis request bit pattern (DISP,VELO, and so on)
    284 AECONFIG(2)  CHAR4 Aerodynamic configuration name
    286 AESYMXY          I Symmetry flag for aerodynamic xy plane
    287 AESYMXZ          I Symmetry flag for aerodynamic xz plane
    288 DISREL           I Displacement relative output flag
    289 VELREL           I Velocity relative output flag
    290 ACCREL           I Acceleration relative output flag
    291 GPEPTSET         I Grid point strain output set (GPSTRAIN)
    292 GPEMEDIA         I Grid point strain output media (GPSTRAIN)
    293 GPEFMT           I Grid point strain output format (GPSTRAIN)
    294 TEMPMAT          I Thermal material set TEMP(MAT).
    295 AECSSSET         I Aerodynamic Control Surface Schedule (CSSCHD)
    296 EKEPTSET         I Element kinetic energy output set (EKE)
    297 EKEMEDIA         I Element kinetic energy media (EKE)
    298 EKEFMT           I Element kinetic energy format (EKE)
    299 EKETHRSH        RS Element kinetic energy threshold (EKE)
    300 EDEPTSET         I Element damping energy output set (EDE)
    301 EDEMEDIA         I Element damping energy media (EDE)
    302 EDEFMT           I Element damping energy format (EDE)
    303 EDETHRSH        RS Element damping energy threshold (EDE)
    304 PANCON           I Panel contributions set (PANCON)
    305 PCMEDIA          I Panel contributions media (PANCON)
    306 PCFMT            I Panel contributions format (PANCON)
    307 PCFORM           I Panel contributions form (PANCON)
    308 PCTOPP           I Panel contributions TOPP (PANCON)
    309 GCTOPG           I Grid contributions TOPG (GRDCON)
    310 PCSOL            I Panel contributions SOLUTION (PANCON)
    311 PCPAN            I Panel contributions PANEL (PANCON)
    312 GCGRID           I Grid contributions GRID (GRDCON)
    313 MODSLF           I Mode selection set (fluid)
    314 EFFMASET         I Modal effective mass output set (MEFFMASS)
    315 EFFMAGID         I Modal effective mass GID (MEFFMASS)
    316 EFFMATHR        RS Modal effective mass fraction threshold (MEFFMASS)
    317 A2GG(2)      CHAR4 Name of direct input (g-set) acoustic coupling matrix (A2GG)
    319 RCRSET           I RCROSS output set
    320 RCRFMT           I RCROSS format
    321 AEUXREF          I AEUXREF
    322 GCHK             I Ground Check Flag (GROUNDCHECK)
    323 GCHKOUT          I Ground Check Output (GROUNDCHECK)
    324 GCHKSET          I Ground Check Set (GROUNDCHECK)
    325 GCHKGID          I Ground Check Gid (GROUNDCHECK)
    326 GCHKTHR         RS Ground Check Thresh (GROUNDCHECK)
    327 GCHKRTHR        RS Ground Check RThresh (GROUNDCHECK)
    328 GCHKDREC         I Ground Check Data recovery (GROUNDCHECK)
    329 ASPCMED          I Output Media Request (AUTOSPC)
    330 ASPCEPS         RS EPS value for fixup (AUTOSPC)
    331 ASPCPRT          I EPS value for printing (AUTOSPC)
    332 ASPCPCH          I Punch Set Id (AUTOSPC)
    333 EXSEGEOM         I External superelement geometry flag (EXTSEOUT)
    334 NA2GG            I Internal set id for A2GG
    335 NK2PP            I Internal set id for K2PP
    336 NM2PP            I Internal set id for M2PP
    337 NB2PP            I Internal set id for B2PP
    338 NK2GG            I Internal set id for K2GG
    339 NM2GG            I Internal set id for M2GG
    340 NB2GG            I Internal set id for B2GG
    341 NP2G             I Internal set id for P2G
    342 GEODSET          I Geometry Check DISP Set identification number (GEOMCHECK)
    343 GEODMXMN         I Geometry Check DISP Max/Min (GEOMCHECK)
    344 GEODOCID         I Geometry Check DISP Max/Min Output Cor. Sys. (GEOMCHECK)
    345 GEODNUMB         I Geometry Check No. of DISP Max/Min Output (GEOMCHECK)
    346 GEOLSET          I Geometry Check OLOAD Set identification number (GEOMCHECK)
    347 GEOLMXMN         I Geometry Check OLOAD Max/Min (GEOMCHECK)
    348 GEOLOCID         I Geometry Check OLOAD Max/Min Output Cor. Sys. (GEOMCHECK)
    349 GEOLNUMB         I Geometry Check No. of OLOAD Max/Min Output (GEOMCHECK)
    350 GEOSSET          I Geometry Check SPCF Set identification number (GEOMCHECK)
    351 GEOSMXMN         I Geometry Check SPCF Max/Min (GEOMCHECK)
    352 GEOSOCID         I Geometry Check SPCF Max/Min Output Cor. Sys. (GEOMCHECK)
    353 GEOSNUMB         I Geometry Check No. of SPCF Max/Min Output (GEOMCHECK)
    354 GEOMSET          I Geometry Check MPCF Set identification number (GEOMCHECK)
    355 GEOMMXMN         I Geometry Check MPCF Max/Min (GEOMCHECK)
    356 GEOMOCID         I Geometry Check MPCF Max/Min Output Cor. Sys. (GEOMCHECK)
    357 GEOMNUMB         I Geometry Check No. of MPCF Max/Min Output (GEOMCHECK)
    358 GEOASET          I Geometry Check ACCE Set identification number (GEOMCHECK)
    359 GEOAMXMN         I Geometry Check ACCE Max/Min (GEOMCHECK)
    360 GEOAOCID         I Geometry Check ACCE Max/Min Output Cor. Sys. (GEOMCHECK)
    361 GEOANUMB         I Geometry Check No. of ACCE Max/Min Output (GEOMCHECK)
    362 GEOVSET          I Geometry Check VELO Set identification number (GEOMCHECK)
    363 GEOVMXMN         I Geometry Check VELO Max/Min (GEOMCHECK)
    364 GEOVOCID         I Geometry Check VELO Max/Min Output Cor. Sys. (GEOMCHECK)
    365 GEOVNUMB         I Geometry Check No. of VELO Max/Min Output (GEOMCHECK)
    366 NTFL             I Internal set id for TFL
    367 BCONTACT         I BCONTACT Set identification number
    368 GPKESET          I Grid point kinetic energy output set (GPKE)
    369 GPKEMEDI         I Grid point kinetic energy media (GPKE)
    370 GPKEFMT          I Grid point kinetic energy format (GPKE)
    371 ELMSUM           I Element Summary Output (ELSUM)
    372 WCHK             I Weight Check Flag (WEIGHTCHECK)
    373 WCHKOUT          I Weight Check Output (WEIGHTCHECK)
    374 WCHKSET          I Weight Check Set identification number (WEIGHTCHECK)
    375 WCHKGID          I Weight Check GID (WEIGHTCHECK)
    376 WCHKCGI          I Weight Check CGI (WEIGHTCHECK)
    377 WCHKWM           I Weight Check Weight/Mass units (WEIGHTCHECK)
    378 EXSEOUT          I External Superelement Output Flag
    379 EXSEMED          I External Superelement Output Media
    380 EXSEUNIT         I External Superelement Output Unit
    381 EXSEASMB         I External Superelement Output ASMBULK Flag
    382 EXSEEXTB         I External Superelement Output EXTBULK Flag
    383 K42GG(2)     CHAR4 Name of direct input (g-set) structural damping matrix K42GG
    385 NK42GG           I Internal set id for K42GG
    386 EXSESTIF         I External Superelement Output STIFFNESS Flag
    387 EXSEMASS         I External Superelement Output MASS Flag
    388 EXSEDAMP         I External Superelement Output DAMPING Flag
    389 EXSEK4DA         I External Superelement Output K4DAMP Flag
    390 EXSELOAD         I External Superelement Output LOADS Flag
    391 EXSESEID         I External Superelement Output SE ID
    392 EXSEDMFX(2)  CHAR4 External Superelement DMIGSFIX String
    394 NSMID            I Non-Structural Mass Set ID
    395 NSELD            I Internal SID for SELOAD
    396 FSELD            I Internal SID for SELOAD scale factor
    397 OP4UNIT          I MBDEXPORT OP4 logical unit number
    398 RPOSTS1          I Random RPOSTS1 parameter
    399 CHECK            I ADAMSMNF/MBDEXPORT CHECK flag
    400 ADMOUT           I ADAMSMNF ADMOUT flag//MBDEXPORT RECVROP2 flag
    401 FLEXBODY         I ADAMSMNF/MBDEXPORT FLEXBODY flag
    402 FLEXONLY         I ADAMSMNF/MBDEXPORT FLEXONLY flag
    403 MINVAR           I ADAMSMNF/MBDEXPORT MINVAR parameter
    404 PSETID           I ADAMSMNF/MBDEXPORT PSETID parameter
    405 OUTGSTRS         I ADAMSMNF/MBDEXPORT OUTGSTRS flag
    406 OUTGSTRN         I ADAMSMNF/MBDEXPORT OUTGSTRN flag
    407 RMSBIT           I Random analysis RMS required bit pattern
    408 MODESCC          I MODES case control existence flag
    409 RMSSF           RS Random RMSSF parameter
    410 UNDEF(3)      None
    413 BCSET            I Contact Set ID
    414 BCRESU           I Contact results output
    415 BCMEDIA          I Contact results media code
    416 BCFMT            I Contact results format code
    417 BCTYPE           I Traction=1, Force=2, Both=3
    418 GKRESU           I Gasket results output
    419 GKMEDIA          I Gasket results media code
    420 GKFMT            I Gasket results format code
    421 PRSSET           I Pressure output set (PRESSURE)
    422 PRSMEDIA         I Pressure output media (PRESSURE)
    423 PRSFMT           I Pressure output format (PRESSURE)
    424 FRFIN            I FRFIN set number
    425 PRSTOTAL         I Pressure output:
      total bit(0)=0,
      scatter bit(0)=0
    426 RSMETHOM         I RSMETHOD parameter
    427 ESETHRSH         I ESE THRESHOLD
    428 MDESET           I Modal energy output set (MODALE)
    429 MDEMEDI          I Modal energy media (MODALE)
    430 MCSOL            I Modal contributions SOLUTION (MODCOM)
    431 MCPAN            I Modal contributions PANELMC (MODCOM)
    432 MDEFMT           I Modal energy output format (MODALE)
    433 ACTLDSET         I Acoustic load set (ALOAD)
    434 MDECMPT          I Modal energy computation set (MODALE)
    435 MDESORT          I Modal energy sort flag (MODALE)
    436 MDETYPE          I Modal energy type flag (MODALE)
    437 MDECALC          I Modal energy calculation flag (MODALE)
    438 RMETSET          I RMETHOD set id
    439 RIGID            I Rigid element type
    440 BOLTPRE          I Bolt preload set
    441 BGSET            I Glue set id
    442 MCTOPF           I Modal contributions TOPF (MODCON)
    443 IPRPU            I RANDOM print/punch option
    444 ADMCHK           I ADMRECVR ADMCHK flag
    445 MODSEL           I Mode selection set (structural)
    446 ADMREC           I ADMRECVR activation flag
    447 ADMFORM          I ADMRECVR ADMFORM parameter
    448 MSRMODE          I ADMRECVR MSRMODE parameter
    449 RGBODY           I ADMRECVR RGBODY flag
    450 MSGLVL           I ADMRECVR MSGLVL parameter
    451 EBDSET           I Element birth/death set
    452 SHELLTHK         I Shell thickness results output flag
    453 STMEDIA          I Shell thickness results media code
    454 STFMT            I Shell thickness results format code
    455 ICTYPE           I Transient IC type
    456 RMXMN            I RMAXMIN flag to indicate presence of card
    457 ROPT             I RMAXMIN print, plot, punch flag
    458 RINP             I RMAXMIN stress, force, displacement flag
    459 RABS             I RMAXMIN maximum, absolute, minimum flag
    460 RAPP             I RMAXMIN approach flag
    461 RMXTRN           I Alternate request of RMXTRN parameter
    462 NPAVG            I Number of maximum peaks to average
    463 RSTAR           RS Start time step for desired interval
    464 RSTOP           RS End time step for desired interval
    465 MODCON           I Modal contribution set
    466 MCMEDIA          I Modal contribution media
    467 MCFMT            I Modal contribution format
    468 MCFORM           I Modal contribution FORM
    469 MCTOPS           I Modal contributions TOPS (MODCON)
    470 PSDD             I SOL200: int. set no. for grids w/ PSDDISP design response
    471 PSDV             I SOL200: int. set no. for grids w/ PSDVELO design response
    472 PSDA             I SOL200: int. set no. for grids w/ PSDACCL design response
    473 ISTAR            I Start subcase id (RMAXMIN)
    474 ISTOP            I End subcase id (RMAXMIN)
    475 FK2PP            I Internal set id for K2PP scale factor
    476 FM2PP            I Internal set id for M2PP scale factor
    477 FB2PP            I Internal set id for B2PP scale factor
    478 FK2GG            I Internal set id for K2GG scale factor
    479 FM2GG            I Internal set id for M2GG scale factor
    480 FB2GG            I Internal set id for B2GG scale factor
    481 FK42GG           I Internal set id for K42GG scale factor
    482 FP2G             I Internal set id for P2G scale factor
    483 FA2GG            I Internal set id for A2GG scale factor
    484 GPRSORT          I Global ply results sorted with global ply ID numbers
    485 EFLOAD1          I External field load orientation
    486 EFLOAD2          I External field coordinate system
    487 BGRESU           I Glue results output
    488 BGMEDIA          I Glue results media code
    489 RANLOOP          I RANDOM loop number; used with ANALYSIS = RANDOM
    490 BGTYPE           I Glue results type
    491 RSVCOMP          I Residual vector component flag
    492 RSVOPTC          I Residual vector component options
    493 RSVSYST          I Residual vector system flag
    494 RSVOPTS          I Residual vector system options
    495 PLSLOC           I Ply strain or stress locations
    496 ELSMOP           I ELSUM output option
    497 ERPSET           I ERP set
    498 ERPSORT          I SORT1/SORT2
    499 ERPMEDIA         I ERP output media
    500 ERPFMT           I ERP output format
    501 ERPSOL           I ERP SOLUTION set
    502 ERPELEM          I ERP element output
    503 ERPCSV           I Unused. Reserved for ERP
    504 ERPCOEFF        RS ERP coefficient
    505 UNDEF(4)      None
    509 ATVFSID          I SID of ATVF
    510 ATVUNIT          I ATVOUT OP2 unit
    511 ATVSETNO         I ATVOUT microphone set identification number
    512 ATVFLAGS         I ATVOUT bits for flags = 1 if ATVOUT specified
    513 ACPANEL          I PANEL in ACPOWER:
      0 for none,
      -1 for all,
      >0 for panel identification number
    514 RMETCMR          I Rotor dynamics CMR method
    515 EFFMAT1         RS Modal effective mass minimum total value in X displacement
    516 EFFMAT2         RS Modal effective mass minimum total value in Y displacement
    517 EFFMAT3         RS Modal effective mass minimum total value in Z displacement
    518 EFFMAMIT         I Modal effective mass maximum number of additional iterations
    519 SEQDEP       CHAR4 Sequence dependency on or off (SEQDEP)
    520 NLCSET           I Set Identification nonlinear control (NLCNTL)
    521 GSTRESET         I Gauss point stress output set (GSTRESS)
    522 GSTMEDIA         I Gauss point stress output media (GSTRESS)
    523 GSTREFMT         I Gauss point stress output format (GSTRESS)
    524 GSTRNSET         I Gauss point total strain output set (GSTRAIN)
    525 GSNMEDIA         I Gauss point total strain output media (GSTRAIN)
    526 GSTRNFMT         I Gauss point total strain output format (GSTRAIN)
    527 ELSTNSET         I Nodal elastic strain on elements, output set (ELSTRN)
    528 ELNMEDIA         I Nodal elastic strain on elements, output media (ELSTRN)
    529 ELSTNFMT         I Nodal elastic strain on elements, output format (ELSTRN)
    530 GELSSET          I Gauss point elastic strains on elements, output set (GELSTRN)
    531 GESMEDIA         I Gauss point elastic strains on elements, output media (GELSTRN)
    532 GELSFMT          I Gauss point elastic strains on elements, output format (GELSTRN)
    533 CRSTSET          I Nodal creep strains on elements, output set (CRSTRN)
    534 CRSMEDIA         I Nodal creep strains on elements, output media (CRSTRN)
    535 CRSTFMT          I Nodal creep strains on elements, output format (CRSTRN)
    536 GCRSSET          I Gauss point creep strains on elements, output set (GCRSTRN)
    537 GCRMEDIA         I Gauss point creep strains on elements, output media (GCRSTRN)
    538 GCRSFMT          I Gauss point creep strains on elements, output format (GCRSTRN)
    539 PLSTSET          I Nodal plastic strains on elements, output set (PLSTRN)
    540 PLSMEDIA         I Nodal plastic strains on elements, output media (PLSTRN)
    541 PLSTFMT          I Nodal plastic strains on elements, output format (PLSTRN)
    542 GPLSSET          I Gauss point plastic strains on elements, output set (GPLSTRN)
    543 GPLMEDIA         I Gauss point plastic strains on elements, output media (GPLSTRN)
    544 GPLSFMT          I Gauss point plastic strains on elements, output format (GPLSTRN)
    545 THSTSET          I Nodal thermal strains on elements, output set (THSTRN)
    546 THSMEDIA         I Nodal thermal strains on elements, output media (THSTRN)
    547 THSTFMT          I Nodal thermal strains on elements, output format (THSTRN)
    548 GTHSSET          I Gauss point thermal strains on elements, output set (GTHSTRN)
    549 GTHMEDIA         I Gauss point thermal strains on elements, output media (GTHSTRN)
    550 GTHSFMT          I Gauss point thermal strains on elements, output format (GTHSTRN)
    551 OTEMPSET         I Temperatures used at solution points, output set (OTEMP)
    552 OTEMEDIA         I Temperatures used at solution points, output media (OTEMP)
    553 OTEMPFMT         I Temperatures used at solution points, output format (OTEMP)
    554 NONCUP           I ADAMSMNF/MBDEXPORT
    555 DTEMPSET         I Time dependent temperature load (DTEMP)
    556 JINSET           I J integral output set (JINTEG)
    557 JINMEDIA         I J integral output media (JINTEG)
    558 ADAPTRESU        I Adaptive Meshing set, output error estimator
    559 ADAPTMEDIA       I Error Estimator media code
    560 ADAPTPYE         I Error Estimator based on ENERGY FORM or STRESS FORM and STEP
    561 INITSSET         I Initial stress/strain:
      INITS=n where n=0 for none,
      n>0 for INITS or INITADD bulk entry SID,
      n<0 for invalid value
    562 OSTNSET          I Set no for initial strain output after subcase 0
    563 OPRESSET         I Pressures used at solution points, output set (OPRESS)
    564 OPRESDIA         I Pressures used at solution points, output media (OPRESS)
    565 BOLTRESMED       I Bolt axial force, shear force, bending moment, and strain output media
    566 CYCLSET          I SOL 401 cyclic symmetry set IF (CYCSET)
    567 OSTNOPT          I Output options for initial strain after subcase 0:
      =1 for element-node;
      =2 for Gauss;
      =3 for both. (For INITSTN/INITSTS)
    568 OSTNMED          I Media for initial strain output after subcase 0; PRINT/PUNCH/PLOT. (For INITSTN/INITSTS)
    569 ACPWRGST         I Acoustic power, GROUP output set (ACPOWER)
    570 ACPWRAST         I Acoustic power, AMLREG output set (ACPOWER)
    571 ACPWRDIA         I Acoustic power, output media (ACPOWER)
    572 ACPWRFMT         I Acoustic power, output format (ACPOWER)
    573 MPINTSET         I Microphone point intensity, output set (ACINTENSITY)
    574 MPINTDIA         I Microphone point intensity, output media (ACINTENSITY)
    575 MPINTFMT         I Microphone point intensity, output format (ACINTENSITY)
    576 OTMFORC          I Output set (OTMFORC)
    577 OTMFORCM         I Output media (OTMFORC)
    578 OTMFORCF         I Output format (OTMFORC)
    579 MPVELSET         I Microphone point velocity, output set (MPVELOCITY)
    580 MPVELDIA         I Microphone point velocity, output media (MPVELOCITY)
    581 MPVELFMT         I Microphone point velocity, output format (MPVELOCITY)
    582 PFRESUSET        I Progressive failure analysis of composites, output set (PFRESULTS)
    583 PFRESUDIA        I Progressive failure analysis of composites, output media (PFRESULTS)
    584 PFRESUFMT        I Progressive failure analysis of composites, output code for damage value/damage status/damage energy (PFRESULTS)
    585 MONVAR           I Maya monitor variable for displacement
    586 CYCFSET          I Forces of cyclic constraint output set (CYCFORCE)
    587 CYCMEDIA         I Forces of cyclic constraint output media (CYCFORCE)
    588 CYCFFMT          I Forces of cyclic constraint output format (CYCFORCE)
    589 BOLTRESULTS      I Bolt axial force, shear force, bending moment, and strain
    590 STVARSET         I State variable values on elements, output set (STATVAR)
    591 STVARMEDIA       I State variable values on elements, output media (STATVAR)
    592 STVARFMT         I State variable values on elements, output format (STATVAR)
    593 CZRESUSET        I Cohesive elements, output set (CZRESULTS)
    594 CZRESUDIA        I Cohesive elements, output media (CZRESULTS)
    595 CZRESUFMT        I Cohesive elements, output code for traction/relative motion/damage value (CZRESULTS)
    596 CKGAPSET         I Gap results, output set (CKGAP)
    597 CKGAPDIA         I Gap results, output media (CKGAP)
    598 CKGAPFMT         I Gap results, output location:
      =1 for grid;
      =2 for Gauss;
      =3 for both (CKGAP)
    599 GRDCON           I Grid contributions set
    600 GCMEDIA          I Grid contributions media
    601 GCFMT            I Grid contributions format
    602 GCFORM           I Grid contributions FORM
    603 GCSOL            I Grid contributions SOLUTION
    604 INITSOFF         I Initial strain offset for balanced initial
      stress/strain: INITS(OFFSET)=n where n=0
      for none, n>0 for INITS or INITADD bulk entry
      SID, n<0 for invalid value
    605 INPWRGST         I Incident acoustic power, GROUP output set (INPOWER)
    606 INPWRFST         I Incident acoustic power, FACES output set (INPOWER)
    607 INPWRDIA         I Incident acoustic power, output media (INPOWER)
    608 INPWRFMT         I Incident acoustic power, output format (INPOWER)
    609 TRPWRGST         I Transmitted acoustic power, GROUP output set (TRPOWER)
    610 TRPWRAST         I Transmitted acoustic power, AMLREG output set (TRPOWER)
    611 TRPWRDIA         I Transmitted acoustic power, output media (TRPOWER)
    612 TRPWRFMT         I Transmitted acoustic power, output format (TRPOWER)
    613 TRLOSFLG         I Acoustic transmission loss, YES/NO flag
      (1=yes, 0=no) (TRLOSS)
    614 TRLOSDIA         I Acoustic transmission loss, output media (TRLOSS)
    615 TRLOSFMT         I Acoustic transmission loss, output format (TRLOSS)
    616 NLARCST          I SOL 401 nonlinear arc-length solution flag set IF (NLARCL)
    617 IMPRFST          I SOL 401 imperfection set flag, SET IF (IMPERF)
    618 MONPNT           I MONPNTn output bit flag(s)
    619 FRFOUT           I Frequency-dependent component output flag (FRFOUT)
    620 FRFOPT           I Frequency-dependent component output options (FRFOUT)
    621 FRFSEID          I SEID for frequency-dependent component output (FRFOUT)
    622 FRFOP2           I Unit for frequency-dependent component output (FRFOUT)
    623 RMSINT           I Random RMSINT parameter
    624 XSEMODAC         I External superelement MODACC parameter
    625 XSEFSCOU         I External superelement FSCOUP parameter
    626 SCSET            I Optimization static subcase set identification number (DESOBJ)
    627 SCFUNC           I Optimization static subcase function option (DESOBJ)
    628 ELAR             I Element add/remove set identification number
    629 ELAROFLG         I Element status output flag:
      1=yes,
      0=no
    (ELAROUT)
    630 DMTRSET           I
      1=yes (default),
      0=no
    631 DMTRMEDIA         I
      bit(1)=1 (default), bit(1)=0 noprint;
      bit(2)=1 punch, bit(2)=0 nopunch (default);
      bit(3)=1 plot
    632 DMTRFMT           I
      0=real/imaginary (default),
      1=magnitude/phase
    633 DMTRTYPE          I Unused
    634 PEAKOUT           I PEAKOUT bulk entry selection
    635 ELAROMDA          I Element status output, output media (ELAROUT)
    636 FLXSLI            I Flexible slider identification number
    637 JCONSET           I Joint constraint set identification number
    638 JRESSET           I Kinematic joints output set (JRESULTS)
    639 JRESMEDIA         I Kinematic joints output media (JRESULTS)
    640 JRESFMT           I Kinematic joints output code:
      1=force,
      2=moment,
      4=position,
      8=rotation,
      16=speed,
      32=rotation speed (JRESULTS)
    641 FLXRSET           I Flexible slider output set (FLXRESULTS)
    642 FLXRMEDIA         I Flexible slider output media (FLXRESULTS)
    643 FLXRFMT           I Flexible slider output code: 64-curvdisp (FLXRESULTS)
    644 ACTEMP            I ACTEMP bulk entry selection
    645 DMTRLSET          I
      1=yes (default),
      0=no
    646 DMTRLSMEDIA       I
       bit(1)=1 print (default), bit(1)=0 noprint;
       bit(2)=1 punch, bit(2)=0 nopunch (default);
       bit(3)=1 plot
    647 DMTRLSFMT         I Unused
    648 ENFUNC            I Optimization entity response function option (DESOBJ)
    649 GPFSOL            I GPFORCE output frequency selection value
    650 CSMSET            I Co-simulation (wetted) region set identification number
    651 DISPSOL           I DISPLACEMENT output frequency selection value
    652 VELOSOL           I VELOCITY output frequency selection value
    653 ACCESOL           I ACCELERATION output frequency selection value
    654 PRESSOL           I PRESSURE output frequency selection value
    655 OPRESSOPT         I Pressures used at solution points (output options):
      0=COUPLED,
      1=FPP,
      2=BOTH
    656 NLCSETG           I Set Identification nonlinear control specified globally (NLCNTL):
      0=none
    657 ACORDCHK      CHAR4 Acoustic max frequency and element order check mode = STOP
    658 UNDEF(542)    None
    LCC LSEM(C)           I Number of symmetry subcase coefficients from item SYMFLG
    The value for LCC is set by word 166
      LCC+1 COEF         RS Symmetry subcase coefficients (SUBSEQ or SYMSEQ)
      Word LCC+1 repeats LSEM times
      LCC+2 SETID         I Set identification number
      LCC+3 SETLEN(C)     I Length of this set
      LCC+4 SETMEM        I Set member identification number
      Word LCC+4 repeats SETLEN times
      Words LCC+2 through LCC+4 repeat NSETS times
      LCC+5 PARA      CHAR4 Hard-coded to "PARA"
      LCC+6 PARLEN(C)     I Length of this parameter value specification
      LCC+7 CHTYPE(C)     I Character type flag: 3 means character, 2 otherwise
      LCC+8 PARAM(2)  CHAR4 Hard-coded to "PARA" and "M "
      LCC+10 PNAME(2) CHAR4 Name of parameter
    PARLEN=8 Length
      LCC+12 INTEGER      I Integer value
    PARLEN=9 Real-double parameter value
      LCC+12 TYPE         I Real type - hard-coded to -4
      LCC+13 REAL        RD Real-double value
    PARLEN=10 Complex-single parameter value
      LCC+12 RTYPE        I Real part type - hard-coded to -2
      LCC+13 REAL        RS Real part value
      LCC+14 ITYPE        I Imaginary part type - hard-coded to -2
      LCC+15 IMAG        RS Imaginary part value
    PARLEN=12 Complex-double parameter value
      LCC+12 RTYPE        I Real part type - hard-coded to -4
      LCC+13 REAL        RD Real part value
      LCC+14 ITYPE        I Imaginary part type - hard-coded to -4
      LCC+15 IMAG        RD Imaginary part value
    End PARLEN
    Words LCC+5 through max repeat until NANQ occurs
    Words LCC+5 through LCC+15 repeat until End of Record
    """
    if nastran_format == 'optistruct':
        raise NotImplementedError(nastran_format)
        #self.show_data(data, types='ifs')
    ints = np.frombuffer(data, dtype=idtype).copy()
    floats = np.frombuffer(data, dtype=fdtype).copy()
    #print(len(ints[15:31]))

    (subcase_id, mpc, spc, load, method_structure, deform,
     temp_load, temp_mat_init, tic,
     nlload_set, nlload_media, nlload_format,
     dload, freq, tfl) = ints[:15]
    (symflag,
     oload_set, oload_media, oload_format,
     disp_set, disp_media, disp_format,
     stress_set, stress_media, stress_format,
     force_set, force_media, force_format,
     accel_set, accel_media, accel_format,
     velocity_set, velocity_media, velocity_format,
     spc_force_set, spc_force_media, spc_force_format,
     tstep) = ints[15:38]

    subcase = Subcase(id=subcase_id)

    title_subtitle_label = data[38*size:134*size]
    title_bytes = title_subtitle_label[:32*size]
    subtitle_bytes = title_subtitle_label[32*size:64*size]
    label_bytes = title_subtitle_label[64*size:]
    #print(f'title    = {title_bytes.strip()!r}')
    #print(f'subtitle = {subtitle_bytes.strip()!r}')
    #print(f'label    = {label_bytes.strip()!r}')
    #--------
    (stpltflag, ax_sym_set, nharmonics, tstrv) = ints[134:138]
    #self.show_data(data[133*size:158*size], types='sdq')
    #print(data[200*size:300*size])
    k2pp_bytes = data[138*size:140*size]
    m2pp_bytes = data[140*size:142*size]
    b2pp_bytes = data[142*size:144*size]
    if size == 8:
        k2pp_bytes = reshape_bytes_block(k2pp_bytes)
        m2pp_bytes = reshape_bytes_block(m2pp_bytes)
        b2pp_bytes = reshape_bytes_block(b2pp_bytes)
    #print(k2pp, m2pp, b2pp)

    (ofreq, sedr, mfluid, cmethod, sdamp) = ints[144:149]
    (sdisp_set, sdisp_media, sdisp_format,
     svelocity_set, svelocity_media, svelocity_format,
     saccel_set, saccel_media, saccel_format,
     nonlinear, partn, cyclic, random, nlparm, fmethod,
     nwords_to_lsem) = ints[150:166]
    if nwords_to_lsem // 4 not in {150, 300}:
        self.log.warning(f'nwords_to_lsem={nwords_to_lsem} nwords_to_lsem//4={nwords_to_lsem//4}')

    #assert nwords_to_lsem == 4 * 150, f'nwords_to_lsem={nwords_to_lsem} size={size} nwords_to_lsem//size={nwords_to_lsem//size}'
    # LCC=1200 @ word 166
    # LSM @ 658?

    #print(ints[144:166].tolist())
    (gpforce_set, gpforce_media, gpforce_format,
     ese_set, ese_media, ese_format,
     aerof_set, aerof_media, aerof_format,
     super_seid, super_load, gust, sefinal) = ints[166:179]

    (svector_set, svector_media, svector_format,
     mpress_set, mpress_media, mpress_format,
     houtput_media1, houtput_media2, houtput_media3,
     houtput_format1, houtput_format2, houtput_format3) = ints[190:202]
    (von_mises_stress, secmd_flag,
     gpstress_set, gpstress_media, gpstress_format,
     strfield_set, strfield_media, strfield_format,
     cload, set2,
     dsa_part, dsa_store, dsa_output,
     strain_set, strain_media, strain_format,
     apressure, trim, modlist, method_fluid,
     elsdcon_set, elsdcon_media, elsdcon_format,
     gpsdcon_set, gpsdcon_media, gpsdcon_format,
     sedv, sere, sers, ) = ints[206:235]

    # 240
    (boutput_set, boutput_media, boutput_format,
     diverg, outrcv, statsub_preload) = ints[235:241]
    adapt, desobj, dessub, subspan, desglb = ints[245:250]
    analysis_bytes = data[250*size:251*size]
    (gpq_stress, gpq_force, gpq_strain,
     suport1, statsub_buckle, bc,
     auxmodel, adact) = ints[251:259]
    (pelement_set, pelement_media, pelement_format,
     vugrid_set, vugrid_media, vugrid_format,
     mpc_force_set, mpc_force_media, mpc_force_format,
     umethod, sdamp_fluid, smethod,
     nlstress_set, nlstress_media, nlstress_format,
     modtrak,
     dsa_form, dsa_expo, dsa_begin, dsa_interval,
     dsa_final, dsa_set) = ints[259:281]
    sort_flag, rand_bit = ints[281:283]
    aeconfig_bytes = data[283*size:285*size]
    aesym_xy, aesym_xz = ints[285:287]

    # ------------
    # 300
    (gpstrain_set, gpstrain_media, gpstrain_format,
     temp_mat, csschd,
     eke_set, eke_media, eke_format,) = ints[290:298]
    eke_thresh = floats[298]
    (ede_set, ede_media, ede_format,) = ints[299:302]
    ede_thresh = floats[302]
    ese_thresh = floats[426]

    meffmass_set, meffmas_grid = ints[313:315]
    meffmass_thresh = floats[315]

    if meffmass_set:
        options = []
        if meffmas_grid:
            options.append(f'GRID={meffmas_grid}')
        if meffmass_thresh:
            options.append(f'THRESH={meffmass_thresh}')
        subcase.add('MEFFMASS', meffmass_set, options, 'STRESS-type')

    autospc_media = ints[328]
    autospc_fixup = floats[329]
    autospc_printing, autospc_punch = ints[330:332]
    nk2pp, nm2pp, nb2pp, nk2gg, nm2gg, nb2gg, np2gg = ints[334:341]
    gpke_set, gpke_media, gpke_format = ints[367:370]

    # we should be able to find PARA and reverse engineer the rest
    # of the table (e.g., MSC 2005/NX 2019)
    #print(data[300*size:10000*size])
    #ddd

    # 3: C:\MSC.Software\simcenter_nastran_2019.2\tpl_post1\cqromidq4.op2
    assert sort_flag in [0, 1, 2, 3], sort_flag  # 3
    #assert rand_bit in [0], rand_bit

    #title = reshape_bytes_block_size(title_bytes, size=size)
    #subtitle = reshape_bytes_block_size(subtitle_bytes, size=size)
    #label = reshape_bytes_block_size(label_bytes, size=size)
    #aeconfig = reshape_bytes_block_size(aeconfig_bytes, size=size)
    #analysis = reshape_bytes_block_size(analysis_bytes, size=size)
    if size == 8:
        title = reshape_bytes_block(title_bytes).decode('latin1').strip()
        subtitle = reshape_bytes_block(subtitle_bytes).decode('latin1').strip()
        label = reshape_bytes_block(label_bytes).decode('latin1').strip()
        aeconfig = reshape_bytes_block(aeconfig_bytes).decode('latin1').strip()
        analysis = reshape_bytes_block(analysis_bytes)# .decode('latin1').strip()
    else:
        title = title_bytes.decode('latin1').strip()
        subtitle = subtitle_bytes.decode('latin1').strip()
        label = label_bytes.decode('latin1').strip()
        aeconfig = aeconfig_bytes.decode('latin1').strip()
        analysis = analysis_bytes # .decode('latin1').strip()
    assert r'\x00' not in str(title), f'{title}'
    assert r'\x00' not in str(subtitle), f'{subtitle}'
    assert r'\x00' not in str(label), f'{label}'
    title = title[:58].rstrip()
    subtitle = subtitle[:55].rstrip()
    label = label[:58].rstrip()
    #print(f'title    = {title!r}')
    #print(f'subtitle = {subtitle!r}')
    #print(f'label    = {label!r}')

    if analysis == b'\x00\x00\x00\x00':
        analysis_str = ''
    elif analysis == b'STAT':
        analysis_str = 'STATICS'
    elif analysis == b'MODE':
        analysis_str = 'MODES'
    elif analysis == b'RAND':
        analysis_str = 'RANDOM'
    elif analysis == b'BUCK':
        analysis_str = 'Buckling'
    elif analysis == b'DYNA':
        analysis_str = 'DYNAMICS'
    elif analysis == b'CYCM':
        analysis_str = 'CYCMODES'
    elif analysis == b'PREL':
        analysis_str = 'PRELOAD'
    elif analysis == b'FOUR':
        analysis_str = 'FOURIER'
    elif analysis == b'TRAN':
        analysis_str = 'TRANSIENT'
    elif analysis == b'SAER':
        analysis_str = 'SAERO'
    elif analysis == b'HEAT':
        analysis_str = 'HEAT'
    elif analysis == b'MTRA':
        analysis_str = 'MTRAN'
    elif analysis == b'MFRE':
        analysis_str = 'MFREQ'
    elif analysis == b'DCEI':
        analysis_str = 'DCEIG'
    elif analysis == b'DFRE':
        analysis_str = 'DFREQ'
    #elif analysis == b'TRAN':
        #analysis_str = 'TRANSIENT'
    #elif analysis == b'TRAN':
        #analysis_str = 'TRANSIENT'
    #elif analysis == b'TRAN':
        #analysis_str = 'TRANSIENT'
    #elif analysis == b'TRAN':
        #analysis_str = 'TRANSIENT'
    else:
        raise NotImplementedError(analysis)

    #sol = 0
    #if analysis_str:
        #sol = 200

    if title:
        subcase.add('TITLE', title, [], 'STRESS-type')
    if subtitle:
        subcase.add('SUBTITLE', subtitle, [], 'STRESS-type')
    if label:
        subcase.add('LABEL', label, [], 'STRESS-type')
    if analysis_str:
        subcase.add('ANALYSIS', analysis_str, [], 'STRESS-type')

    #print('%r' % title)
    #print('%r' % subtitle)
    #print('%r' % label)
    #print('%r' % aeconfig)

    # same up to ~383
    #qqq
    if apressure == -1:
        apressure = 'ALL'
    if sedr == -1:
        sedr = 'ALL'

    if desobj < 0:
        desobj = 0

    basic = {
        #'SOL' : sol,
        'MPC' : mpc,
        'SPC' : spc,
        'LOAD' : load,
        'DEFORM' : deform,
        'TEMP(LOAD)' : temp_load,
        'TEMP(INIT)' : temp_mat_init,
        'TEMP(MAT)' : temp_mat,
        'CSSCHD' : csschd,
        'IC' : tic,
        'DLOAD': dload,
        'FREQ' : freq,
        'TFL' : tfl,
        'TSTEP' : tstep,
        'SEDR' : sedr,
        'MFLUID' : mfluid,
        'METHOD' : method_structure,
        'CMETHOD' : cmethod,
        'TRIM' : trim,
        'RANDOM' : random,
        'STATSUB(PRELOAD)' : statsub_preload,
        'STATSUB(BUCKLE)' : statsub_buckle,
        'DESOBJ' : desobj,
        'DESSUB' : dessub,
        'DRSPAN' : subspan,
        'DESGLB' : desglb,
        'SUPORT1' : suport1,
        'BC' : bc,
        'SDAMP(FLUID)' : sdamp_fluid,
        'SMETHOD' : smethod,
        'AESYMXY' : aesym_xy,
        'AESYMXZ' : aesym_xz,
        'AXISYMMETRIC' : ax_sym_set,
        'HARMONICS' : nharmonics,
        'APRESSURE' : apressure,
        'DIVERG' : diverg,
        'OUTRCV' : outrcv,
        'ADAPT' : adapt,
    }
    stress_options = {
        'NLLOAD': (nlload_set, nlload_media, nlload_format),
        'OLOAD' : (oload_set, oload_media, oload_format),

        'DISP' : (disp_set, disp_media, disp_format),
        'VELOCITY' : (velocity_set, velocity_media, velocity_format),
        'ACCEL' : (accel_set, accel_media, accel_format),

        'FORCE' : (force_set, force_media, force_format),
        'STRESS' : (stress_set, stress_media, stress_format),
        'STRAIN' : (strain_set, strain_media, strain_format),
        'SPCFORCE' : (spc_force_set, spc_force_media, spc_force_format),
        'GPFORCE' : (gpforce_set, gpforce_media, gpforce_format),
        'ESE' : (ese_set, ese_media, ese_format),
        'EKE' : (eke_set, eke_media, eke_format),
        'EDE' : (ede_set, ede_media, ede_format),
        'GPKE' : (gpke_set, gpke_media, gpke_format),

        'ELSDCON' : (elsdcon_set, elsdcon_media, elsdcon_format),
        'GPSDCON' : (gpsdcon_set, gpsdcon_media, gpsdcon_format),
        'GPSTRESS' :(gpstress_set, gpstress_media, gpstress_format),
        'STRFIELD' :(strfield_set, strfield_media, strfield_format),
        'GPSTRAIN' : (gpstrain_set, gpstrain_media, gpstrain_format),

        'SDISP' : (sdisp_set, sdisp_media, sdisp_format),
        'SVELO' : (svelocity_set, svelocity_media, svelocity_format),
        'SACCEL' : (saccel_set, saccel_media, saccel_format),
        'SVECTOR' : (svector_set, svector_media, svector_format),
        'DATAREC' : (pelement_set, pelement_media, pelement_format),
        'VUGRID' : (vugrid_set, vugrid_media, vugrid_format),
        'MPCFORCE' : (mpc_force_set, mpc_force_media, mpc_force_format),
        'NLSTRESS' : (nlstress_set, nlstress_media, nlstress_format),
        'BOUTPUT' : (boutput_set, boutput_media, boutput_format),
    }
    names_namesi = [
        ('K2PP', k2pp_bytes),
        ('M2PP', m2pp_bytes),
        ('B2PP', b2pp_bytes),
     ]
    for key, name_bytes in names_namesi:
        if name_bytes == b'\x00\x00\x00\x00\x00\x00\x00\x00':
            continue
        name = name_bytes.decode('latin1').strip()
        subcase.add(key, name, [], 'STRESS-type')

    for key, value in basic.items():
        if value == 0:
            continue
        #print(key, value)
        subcase.add(key, value, [], 'STRESS-type')

    for key, value in stress_options.items():
        (value_set, value_media, value_format) = value
        if value_set == 0:
            continue
        if value_set == -1:
            value_set = 'ALL'

        #print('***', key, value_set, value_media, value_format)

        if value_media == 0:
            media = []
        elif value_media == 1:
            media = ['PRINT']
        elif value_media == 2:
            media = ['PLOT']
        elif value_media == 3:
            media = ['PLOT', 'PRINT']
        elif value_media == 4:
            media = ['PUNCH']
        elif value_media == 5:
            media = ['PRINT', 'PUNCH']
        elif value_media == 6:
            media = ['PLOT', 'PUNCH']
        elif value_media == 7:
            media = ['PLOT', 'PRINT', 'PUNCH']
        elif value_media == 17:
            media = ['PRINT', 'AVERAGE']
        else:
            print(subcase)
            raise NotImplementedError((key, value_media))

        default_sort = []
        if sort_flag == 0:
            default_sort = []
        elif sort_flag == 1:
            default_sort = ['SORT1']
        elif sort_flag == 2:
            default_sort = ['SORT2']
        #else:
            #raise NotImplementedError(default_sort)

        if value_format == 0:
            sort_method = default_sort
        elif value_format in [1]:
            sort_method = ['SORT1']
        elif value_format in [-1]:
            sort_method = ['SORT2']
        elif value_format in [2]:
            sort_method = ['IMAG'] + default_sort
        elif value_format in [-2]:
            sort_method = ['SORT2', 'IMAG']
        elif value_format == 3:
            sort_method = ['PHASE'] + default_sort
        elif value_format == -3:
            sort_method = ['SORT2', 'PHASE']
        #elif value_format == -1:
            #sort_method = ['SORT2']
        else:
            print(subcase)
            raise NotImplementedError((key, value_format))
        options = media + sort_method
        if key == 'ESE' and ese_thresh:
            options.append(f'THRESH={ese_thresh}')
        elif key == 'EKE' and eke_thresh:
            options.append(f'THRESH={eke_thresh}')
        elif key == 'EDE' and ede_thresh:
            options.append(f'THRESH={ede_thresh}')

        subcase.add(key, value_set, options, 'STRESS-type')
    return subcase
