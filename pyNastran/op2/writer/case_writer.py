from __future__ import annotations
from copy import deepcopy
from struct import pack, Struct
from collections import defaultdict
from typing import TYPE_CHECKING

from pyNastran.op2.errors import SixtyFourBitError
from .geom1_writer import write_geom_header, close_geom_table, init_table
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.subcase import Subcase

def write_casecc_header(table_name: bytes, op2_file, op2_ascii, endian: bytes=b'<'):
    op2_ascii.write('----------\n')
    data = init_table(table_name)
    op2_file.write(pack('4i 8s i 3i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        4, 7, 4,
        #28, 1, 2, 3, 4, 5, 6, 7, 28,

        28,
        101, 1, 0, 1030, 0, 0, 0,
        28,
    ]
    #struct_3i = Struct(endian + b'3i')
    op2_file.write(pack('3i 9i', *data))
    op2_ascii.write(str(data) + '\n')

    #-------------------------------------
    data = [
        4, -2, 4,
        4, 1, 4,
        4, 0, 4]
    op2_file.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

    data = [
        #4, 0, 4,
        4, 2, 4,
        8, b'XCASECC ', 8,
    ]
    op2_file.write(pack('3i i8si', *data))
    op2_ascii.write(str(data) + '\n')
    #data = [8, 1, 2, 8]
    #op2_file.write(pack('4i', *data))
    #-------------------------------------

    data = [
        4, -3, 4,
        4, 1, 4,
        4, 0, 4]
    op2_file.write(pack('9i', *data))
    op2_ascii.write(str(data) + '\n')

def write_casecc(op2_file, op2_ascii, obj, endian: bytes=b'<',
                 nastran_format: str='nx'):
    """writes the CASECC table"""
    write_casecc_header(b'CASECC', op2_file, op2_ascii, endian=endian)

    itable = -3
    subcases = obj.subcases
    if nastran_format == 'msc':
        for subcase_id, subcase in sorted(subcases.items()):
            print(subcase)
            write_msc_casecc(subcase_id, subcase, obj)
    else:
        asdf
    close_geom_table(op2_file, op2_ascii, itable)

def _get_int(key: str, subcase: Subcase) -> int:
    value = 0
    if key in subcase:
        value = subcase[key][0]
        assert isinstance(value, int), type(value)
    return value

def _get_str(key: str, subcase: Subcase, nbytes: int) -> bytes:
    assert nbytes > 0, nbytes
    if key in subcase:
        value_str = subcase[key][0]
        assert isinstance(value_str, str), value_str
        fmt = '%-%%ss' % nbytes
        value_str2 =  fmt % value_str
        value_bytes = value_str2.encode('ascii')
        assert len(v)
    else:
        value_bytes = b' ' * nbytes
    return value_bytes

def _get_stress(key: str, subcase: Subcase) -> int:
    value = 0
    media = 0
    fmt = 0
    von_mises = 1
    is_fmt = False
    if key in subcase:
        value, options_ = subcase[key]
        options = deepcopy(options_)
        if 'SORT1' in options:
            options.remove('SORT1')
            assert is_fmt is False
            fmt = 1
            is_fmt = True
        if 'SORT2' in options:
            options.remove('SORT2')
            assert is_fmt is False
            fmt = 2

        if 'PLOT' in options:
            options.remove('PLOT')
            media += 1
        if 'PRINT' in options:
            options.remove('PRINT')
            media += 2
        if 'PUNCH' in options:
            options.remove('PUNCH')
            media += 4
        assert len(options) == 0, options
    return value, media, fmt, von_mises

def _get_set_media_load(key: str, subcase: Subcase) -> int:
    #if media in (1,    3,    5,    7):
        #options.append('PLOT')
    #if media in (   2, 3,       6, 7):
        #options.append('PRINT')
    #if media in (         4, 5, 6, 7):
        #options.append('PUNCH')
    value = 0
    media = 0
    fmt = 0
    is_fmt = False
    if key in subcase:
        value, options_ = subcase[key]
        options = deepcopy(options_)
        if 'SORT1' in options:
            options.remove('SORT1')
            assert is_fmt is False
            fmt = 1
            is_fmt = True
        if 'SORT2' in options:
            options.remove('SORT2')
            assert is_fmt is False
            fmt = 2

        if 'PLOT' in options:
            options.remove('PLOT')
            media += 1
        if 'PRINT' in options:
            options.remove('PRINT')
            media += 2
        if 'PUNCH' in options:
            options.remove('PUNCH')
            media += 4
        assert len(options) == 0, options
    return value, media, fmt

def write_msc_casecc(subcase_id: int, subcase: Subcase, model: BDF):
    nsub = 0
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
    """
    mpc_id = _get_int('MPC', subcase)
    spc_id = _get_int('SPC', subcase)
    load_id = _get_int('LOAD', subcase)
    deform_id = _get_int('DEFORM', subcase)
    temp_load_id = _get_int('TEMP(LOAD)', subcase)
    temp_mat_id = _get_int('TEMP(MAT)', subcase)
    ic_id = _get_int('IC', subcase)
    ree_set = _get_int('METHOD', subcase)

    nlload_set, nlload_media, nlload_fmt = _get_set_media_load('NLLOAD', subcase)

    data = [subcase_id, mpc_id, spc_id, load_id, ree_set, deform_id, temp_load_id, temp_mat_id, ic_id,
            nlload_set, nlload_media, nlload_fmt]
    """
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
    """
    symflag = 0
    oload_set, oload_media, oload_format = _get_set_media_load('OLOAD', subcase)
    disp_set, disp_media, disp_format = _get_set_media_load('DISPLACEMENT', subcase)
    stress_set, stress_media, stress_set = _get_set_media_load('STRESS', subcase)
    force_set, force_media, force_set = _get_set_media_load('FORCE', subcase)

    data += [symflag,
             oload_set, oload_media, oload_format,
             disp_set, disp_media, disp_format,
             stress_set, stress_media, stress_set,
             force_set, force_media, force_set]

    """
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
    """
    n32 = 7
    nsub -= (3 + n32 * 3)
    accel_set, accel_media, accel_fmt = _get_set_media_load('ACCELERATION', subcase)
    velo_set, velo_media, velo_fmt = _get_set_media_load('VELOCITY', subcase)
    spc_force_set, spc_force_media, spc_force_fmt = _get_set_media_load('SPCFORCE', subcase)
    #oload_set, oload_media, oload_format = _get_set_media_load('MPCFORCE', subcase)
    sdisp_set, sdisp_media, sdisp_format = _get_set_media_load('SDISPLACEMENT', subcase)
    svelo_set, svelo_media, svelo_fmt = _get_set_media_load('SVELOCITY', subcase)
    saccel_set, saccel_media, saccel_fmt = _get_set_media_load('SACCELERATION', subcase)

    tstep_id = _get_int('TSTEP', subcase)
    title = _get_str('TITLE', subcase, nbytes=128)
    subtitle = _get_str('SUBTITLE', subcase, nbytes=128)
    label = _get_str('LABEL', subcase, nbytes=128)
    model_plot_flag = 0
    axisymmetric_flag = 0
    nharmonics = 0
    needs_definition = 0
    k2pp_name = _get_str('K2PP', subcase, nbytes=8)
    m2pp_name = _get_str('M2PP', subcase, nbytes=8)
    b2pp_name = _get_str('B2PP', subcase, nbytes=8)
    ofreq_otime_id = 0
    sedr = _get_int('SEDR', subcase)
    mfluid_id = _get_int('SEDR', subcase)
    cmethod_id = _get_int('SEDR', subcase)
    sdamp_id = _get_int('SDAMP', subcase)

    data += [accel_set, accel_media, accel_fmt,
             velo_set, velo_media, velo_fmt,
             spc_force_set, spc_force_media, spc_force_fmt,
             tstep_id, title, subtitle, label, model_plot_flag, axisymmetric_flag, nharmonics, needs_definition,
             k2pp_name, m2pp_name, b2pp_name, ofreq_otime_id, sedr, mfluid_id, cmethod_id, sdamp_id,
             sdisp_set, sdisp_media, sdisp_format,
             svelo_set, svelo_media, svelo_fmt,
             saccel_set, saccel_media, saccel_fmt,
             ]
    #assert len(data) == 159 - nsub, len(data)
    """
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
    """
    nonlinear_id = _get_int('NONLINEAR', subcase)
    partn = 0
    dsym_id = 0
    random_id = _get_int('RANDOM', subcase)
    nlparm_id = _get_int('NLPARM', subcase)
    fmethod_id = _get_int('FMETHOD', subcase)
    nwords_to_lsem = -999999
    gpforce_set, gpforce_media, gpforce_format = _get_set_media_load('GPFORCE', subcase)
    ese_set, ese_media, ese_format = _get_set_media_load('ESE', subcase)
    aerof_set, aerof_media, aerof_format = _get_set_media_load('AEROF', subcase)
    data += [nonlinear_id, partn, dsym_id, random_id, nlparm_id, fmethod_id,
             nwords_to_lsem,
             gpforce_set, gpforce_media, gpforce_format,
             ese_set, ese_media, ese_format,
             aerof_set, aerof_media, aerof_format]
    #assert len(data) == 175 - nsub, len(data)
    """
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
    """
    nsub -= 3
    seid = _get_int('SUPER', subcase)
    super_id = _get_int('SUPER', subcase)
    gust_id = _get_int('GUST', subcase)
    sefinal_id = _get_int('SEFINAL', subcase)
    semg = _get_int('SEMG', subcase)
    sekr = _get_int('SEKR', subcase)
    selg = _get_int('SELG', subcase)
    selr = _get_int('SELR', subcase)
    seexclude = _get_int('SEEXCLUDE', subcase)
    k2gg = _get_str('K2GG', subcase, nbytes=8)
    m2gg = _get_str('M2GG', subcase, nbytes=8)
    b2gg = _get_str('B2GG', subcase, nbytes=8)
    data += [seid, super_id, gust_id, sefinal_id, semg, sekr, selg, selr, seexclude,
             k2gg, m2gg, b2gg]
    #assert len(data) == 189 - nsub, len(data)

    """
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
    """
    svector_set, svector_media, svector_fmt = _get_set_media_load('SVECTOR', subcase)
    mpress_set, mpress_media, mpress_fmt = _get_set_media_load('MPRES', subcase)
    p2g = _get_str('P2GG', subcase, nbytes=8)
    loadset_id = _get_int('LOADSET', subcase)
    #semg = _get_int('SEMG', subcase)
    semr = _get_int('SEMR', subcase)
    stress_set, stress_media, stress_fmt, von_mises = _get_stress('STRESS', subcase)
    houtput = 0
    noutput = 0
    se_cmd_flag = 0
    data += [
        svector_set, svector_media, svector_fmt,
        mpress_set, mpress_media, mpress_fmt,
        houtput, houtput, houtput,
        noutput, noutput, noutput,
        p2g, loadset_id, semr, von_mises, se_cmd_flag
        #svector_set, svector_media, svector_fmt,
    ]
    #assert len(data) == 208, len(data)
    """
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
    """
    gpstress_set, gpstress_media, gpstress_fmt = _get_set_media_load('GPSTRESS', subcase)
    strfield_set, strfield_media, strfield_fmt = _get_set_media_load('STRFIELD', subcase)
    #sensity_set, sensity_media, sensity_fmt = _get_set_media_load('SENSITY', subcase)
    strain_set, strain_media, strain_fmt = _get_set_media_load('STRAIN', subcase)
    cload_id = _get_int('CLOAD', subcase)
    set2_id = 0
    dsa_print_id = 0
    dsa_store_id = 0
    dsa_output_id = 0
    data += [
        gpstress_set, gpstress_media, gpstress_fmt,
        strfield_set, strfield_media, strfield_fmt,
        cload_id, set2_id,
        dsa_print_id, dsa_store_id, dsa_output_id,
        strain_set, strain_media, strain_fmt,
    ]
    #assert len(data) == 222, len(data)
    """
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
    """
    apress = _get_int('APRESSURE', subcase)
    trim_id = _get_int('TRIM', subcase)
    omodes = _get_int('OMODES', subcase)
    method_fluid = _get_int('METHOD(FLUID)', subcase)
    elsdcon_set, elsdcon_media, elsdcon_fmt = _get_set_media_load('ELSDCON', subcase)
    gpsdcon_set, gpsdcon_media, gpsdcon_fmt = _get_set_media_load('GPSDCON', subcase)
    data += [
        apress, trim_id, omodes, method_fluid,
        elsdcon_set, elsdcon_media, elsdcon_fmt,
        gpsdcon_set, gpsdcon_media, gpsdcon_fmt,
    ]
    """
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
    """
    sedv = _get_int('SEDV', subcase)
    sere = _get_int('SERE', subcase)
    sers = _get_int('SERS', subcase)
    boutput_set, boutput_media, boutput_fmt = _get_set_media_load('BOUTPUT', subcase)
    diverg_id = _get_int('DIVERG', subcase)
    outrcv = 0
    statsub_preload = _get_int('STATSUB(PRELOAD)', subcase)
    modes_select_structure = 0
    mode_select_fluid = 0

    sol = 0
    if model.sol is not None:
        sol = model.sol

    sol_method = 0
    if model.sol_method is not None:
        sol_method = model.sol_method

    data += [
        sedv, sere, sers,
        boutput_set, boutput_media, boutput_fmt,
        diverg_id, outrcv,
        statsub_preload, modes_select_structure, mode_select_fluid,
        sol, sol_method,
    ]
    """
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
    """
    adapt = _get_int('ADAPT', subcase)
    desobj = _get_int('DESOBJ', subcase)
    dessub = _get_int('DESSUB', subcase)
    subspan = _get_int('DRSPAN', subcase)
    desglb = _get_int('DESGLB', subcase)
    analysis = _get_str('ANALYSIS', subcase, nbytes=4)
    stress_corner, force_corner, strain_corner = 0, 0, 0
    suport1 = _get_int('SUPORT1', subcase)
    statsub_buckle = _get_int('STATSUB(BUCKLE)', subcase)
    bc = _get_int('BC', subcase)

    data += [
        adapt, desobj, dessub, subspan, desglb, analysis,
        stress_corner, force_corner, strain_corner,
        suport1, statsub_buckle, bc,
    ]
    """
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
    """
    auxmodel = _get_int('AUXMODEL', subcase)
    adact = _get_int('ADACT', subcase)
    data_set, data_media, data_fmt = _get_set_media_load('DATAREC', subcase)
    vu_set, vu_media, vu_fmt = _get_set_media_load('VUGRID', subcase)
    mpc_set, mpc_media, mpc_fmt = _get_set_media_load('MPCFORCE', subcase)
    nlstress_set, nlstress_media, nlstress_fmt = _get_set_media_load('NLSTRESS', subcase)
    umethod = 0
    sdamp_fluid = _get_int('SDAMP(FLUID)', subcase)
    smethod = _get_int('SMETHOD', subcase)
    data += [
        auxmodel, adact,
        data_set, data_media, data_fmt,
        vu_set, vu_media, vu_fmt,
        mpc_set, mpc_media, mpc_fmt,
        umethod, sdamp_fluid, smethod,
        nlstress_set, nlstress_media, nlstress_fmt,
    ]
    """
    275 MODTRKID         I Mode tracking control parameter set (MODTRAK)
    276 DSAFORM          I Design sensitivity output format: 1=yes,2=no (DSAPRT)
    277 DSAEXPO          I Design sensitivity output export: 1=no,2=yes (DSAPRT)
    278 DSABEGIN         I Design sensitivity output start iteration (DSAPRT)
    279 DSAINTVL         I Design sensitivity output interval (DSAPRT)
    280 DSAFINAL         I Design sensitivity output final iteration (DSAPRT)
    281 DSASETID         I Design sensitivity output set (DSAPRT)
    282 SORTFLG          I Overall SORT1/SORT2 flag: 1 means SORT1 and 2 means SORT2.
    283 RANDBIT          I Random analysis request bit pattern (DISP,VELO, and so on)
    """
    modtrack_id = _get_int('MODTRAK', subcase)
    dsa_form, dsa_expo, dsa_begain, dsa_interval, dsa_final, dsa_set = 0, 0, 0, 0, 0, 0
    sort_flag = 1
    randbit = 0
    data += [
        modtrack_id, dsa_form, dsa_expo, dsa_begain, dsa_interval, dsa_final, dsa_set,
        sort_flag, randbit,
    ]
    """
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
    """
    aeconfig = b'AECONFIG'
    ae_sym_xy = _get_int('AESYMXY', subcase)
    ae_sym_xz = _get_int('AESYMXZ', subcase)
    disp_rel, velo_rel, accel_rel = 0, 0, 0
    temp_mat_id2 = _get_int('TEMP(MAT)', subcase)
    csschd_id = _get_int('CSSCHD', subcase)
    gpstrain_set, gpstrain_media, gpstrain_fmt = _get_set_media_load('GPSTRAIN', subcase)
    eke_set, eke_media, eke_fmt = _get_set_media_load('EKE', subcase)
    ede_set, ede_media, ede_fmt = _get_set_media_load('EDE', subcase)
    eke_threshold = ede_threshold = 0.0
    data += [
        aeconfig, ae_sym_xy, ae_sym_xz,
        disp_rel, velo_rel, accel_rel,
        gpstrain_set, gpstrain_media, gpstrain_fmt,
        temp_mat_id2, csschd_id,
        eke_set, eke_media, eke_fmt, eke_threshold,
        ede_set, ede_media, ede_fmt, ede_threshold,
    ]
    """
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
    """
    pancon_set, pancon_media, pancon_fmt = _get_set_media_load('PANCON', subcase)
    gridcon, mode_select = 0, 0
    data += [
        pancon_set, pancon_media, pancon_fmt, pancon_form, pancon_topp, pancon_topg, pancon_solution, pancon_panel,
        gridcon, mode_select]
    """
    314 EFFMASET         I Modal effective mass output set (MEFFMASS)
    315 EFFMAGID         I Modal effective mass GID (MEFFMASS)
    316 EFFMATHR        RS Modal effective mass fraction threshold (MEFFMASS)
    317 A2GG(2)      CHAR4 Name of direct input (g-set) acoustic coupling matrix (A2GG)
    319 RCRSET           I RCROSS output set
    320 RCRFMT           I RCROSS format
    321 AEUXREF          I AEUXREF
    """
    meffmass_set, meffmass_node, meffmasss_threshold = _get_set_media_load('MEFFMASS', subcase)
    a2gg_name = _get_str('A2GG', subcase, nbytes=8)
    rcross_set, rcross_media, rcross_fmt = _get_set_media_load('RCROSS', subcase)
    aeuxref = 0
    data += [
        meffmass_set, meffmass_node, meffmasss_threshold,
        a2gg_name,
        rcross_set, rcross_fmt,
        aeuxref]
    """
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
    """
    data = []
    """
    333 EXSEGEOM         I External superelement geometry flag (EXTSEOUT)
    334 NA2GG            I Internal set id for A2GG
    335 NK2PP            I Internal set id for K2PP
    336 NM2PP            I Internal set id for M2PP
    337 NB2PP            I Internal set id for B2PP
    338 NK2GG            I Internal set id for K2GG
    339 NM2GG            I Internal set id for M2GG
    340 NB2GG            I Internal set id for B2GG
    341 NP2G             I Internal set id for P2G
    """
    """
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
    """
    """
    366 NTFL             I Internal set id for TFL
    367 BCONTACT         I BCONTACT Set identification number
    368 GPKESET          I Grid point kinetic energy output set (GPKE)
    369 GPKEMEDI         I Grid point kinetic energy media (GPKE)
    370 GPKEFMT          I Grid point kinetic energy format (GPKE)
    371 ELMSUM           I Element Summary Output (ELSUM)
    """
    """
    372 WCHK             I Weight Check Flag (WEIGHTCHECK)
    373 WCHKOUT          I Weight Check Output (WEIGHTCHECK)
    374 WCHKSET          I Weight Check Set identification number (WEIGHTCHECK)
    375 WCHKGID          I Weight Check GID (WEIGHTCHECK)
    376 WCHKCGI          I Weight Check CGI (WEIGHTCHECK)
    377 WCHKWM           I Weight Check Weight/Mass units (WEIGHTCHECK)
    """
    """
    378 EXSEOUT          I External Superelement Output Flag
    379 EXSEMED          I External Superelement Output Media
    380 EXSEUNIT         I External Superelement Output Unit
    381 EXSEASMB         I External Superelement Output ASMBULK Flag
    382 EXSEEXTB         I External Superelement Output EXTBULK Flag
    """
    """
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
    """
    """
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
    """
    """
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
    """
    data += []
    """
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
    asdf
    assert -999999 not in data
