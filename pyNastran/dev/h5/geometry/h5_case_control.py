from __future__ import annotations
from copy import deepcopy
from typing import TYPE_CHECKING
from pyNastran.bdf.case_control_deck import CaseControlDeck
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF
    from pyNastran.bdf.subcase import Subcase
    from cpylog import SimpleLogger

def load_parameters(hdf5_file, geom_model: BDF) -> None:
    """reads PARAMs"""
    pvt_str = 'NASTRAN/INPUT/PARAMETER/PVT'
    parameters = hdf5_file.get(pvt_str)
    for name in list(parameters):
        groupi = parameters[name]
        assert len(groupi.dtype.names) == 3, groupi.dtype.names
        NAME = groupi['NAME']
        DOMAIN_ID = groupi['DOMAIN_ID']
        if name in ('INT', 'DOUBLE'):
            #('NAME', 'VALUE', 'DOMAIN_ID')
            VALUE = groupi['VALUE']
            for key_bytes, value in zip(NAME, VALUE):
                key = key_bytes.strip().decode('latin1')
                obj = geom_model.add_param(key, value, comment='')
                str(obj)
        elif name == 'CHAR':
            #('NAME', 'VALUE', 'DOMAIN_ID')
            VALUE = groupi['VALUE']
            for key_bytes, value_bytes in zip(NAME, VALUE):
                key = key_bytes.strip().decode('latin1')
                value = value_bytes.strip().decode('latin1')
                obj = geom_model.add_param(key, value, comment='')
                str(obj)
        else:
            raise RuntimeError(name)


def load_case_control(geom_model: BDF, hdf5_case_control):
    """reads the case control deck"""
    names = hdf5_case_control.dtype.names
    names = set(list(names))
    names.remove('SID')
    lines = []
    cc = CaseControlDeck(lines, log=None)
    log = cc.log
    sids = hdf5_case_control['SID']
    log.info(f'subcases = {sids}')
    for sid in sids:
        subcase = cc.create_new_subcase(sid)
        cc.subcases[sid] = subcase

    integer_mapped_names = {
        'NOHARMON': ('HARMONICS', 1),

        'SPCSET': ('SPC', 0),
        'MPCSET': ('MPC', 0),
        'TIC': ('IC', 0),
        'ESLSET': ('LOAD', 0),
        'DYMLDSET': ('DLOAD', 0),
        'TFSET': ('TFL', 0),
        'FLUTTER': ('FLUTTER', 0),
        'RANDOM': ('RANDOM', 0),
        'GUST': ('GUST', 0),
        'LOADSET': ('LOADSET', 0),
        'CLOAD': ('CLOAD', 0),
        'APRESS': ('APRESS', 0),
        'TRIM': ('TRIM', 0),

        'FLDBNDY': ('MFLUID', 0),
        'DIVERG': ('DIVERG', 0),
        'DESOBJ': ('DESOBJ', 0),
        'DESSUB': ('DESSUB', 0),
        'BCID': ('BC', 0),
        'REESET': ('METHOD(STRUCTURE)', 0),
        'ELDSET': ('DEFORM', 0),
        'THLDSET': ('TEMP(LOAD)', 0),
        'THMATSET': ('TEMP(MAT)', 0),  # or TEMP(INIT) ???
        'FEQRESET': ('FREQUENCY', 0),
        'TSTEPTRN': ('TSTEP', 0),
        'AXSYMSET': ('AXISYMMETRIC', 0),
        'DAMPTBL': ('SDAMP(STRUCT)', 0),
        'NONPARAM': ('NLPARM', 0),
        'LOADSET': ('LOADSET', 0),
        'MODLIST': ('OMODES', 0),
        'REESETF': ('METHOD(FLUID)', 0),
        'NMODES': ('NMODES', 0),
        'ADAPT' : ('ADAPT', 0),
        'SUBSPAN' : ('SUBSPAN', 0),
        'DESGLB' : ('DESGLB', 0),
        'SUPORT1' : ('SUPORT1', 0),
        'STATSUBB' : ('STATSUB(BUCKLE)', 0),
        'STATSUBP' : ('STATSUB(PRELOAD)', 0),
        'AUXMODEL' : ('AUXMODEL', 0),
        'NSMID': ('NSM', 0),
        'NONLINLD': ('NONLINEAR', 0),
        'CYCLIC': ('DSYM', 0),
        'DYNRED': ('DYNRED', 0),
        'CEESET': ('CMETHOD', 0),
        'AECSSSET': ('CSSCHD', 0),
        'TEMPMAT': ('TEMP(MAT)', 0),
        'RGYRO': ('RGYRO', 0),

        # superelements
        'SEMG' : ('SEMG', -1),
        'SEKR' : ('SEKR', -1),
        'SELG' : ('SELG', -1),
        'SELR' : ('SELR', -1),
        #'SEID' : ('SEID', -1),
        'SEMR' : ('SEMR', -1),
        'SEDV' : ('SEDV', -1),
        'SERE' : ('SERE', -1),
        'SEEXCLUDE' : ('SEEX', 0),
    }

    stress_mapped_names = [
        # media = plot; print; punch
        # format = real; real/imag; magnitude/phase
        # rout - random output
        # RANDBIT: Random analysis request bit pattern (DISP,VELO,etc.)
        ('DISPLACEMENT', 'DPLPTSET', 'DPLMEDIA', 'DPLFMT'), # ROUTDISP
        ('OLOAD', 'LDSPTSET', 'LDSMEDIA', 'LDSFMT'), # ROUTLOAD
        ('VELOCITY', 'VELPTSET', 'VELMEDIA', 'VELFMT'), # ROUTVELO
        ('ACCEL', 'ACCPTSET', 'ACCMEDIA', 'ACCFMT'), # ROUTACCE
        ('MPCFORCE', 'MPCFSET', 'MPCMEDIA', 'MPCFFMT'), # ROUTMSCF
        ('SPCFORCE', 'FOCPTSET', 'FOCMEDIA', 'FOCFMT'), # ROUTSPCF

        ('FORCE', 'FCEPTSET', 'FCEMEDIA', 'FCEFMT'), # ROUTFORC
        ('STRESS', 'STSPTSET', 'STSMEDIA', 'STSFMT'), # ROUTSTRS, VONMISES-done
        ('STRAIN', 'STNSET', 'STNMEDIA', 'STNFMT'), # ROUTSTRN
        #'GPQSTRS', # CQUAD4 corner stress output - done
        #'GPQFORC', # CQUAD4 corner force output - done
        #'GPQSTRN', # CQUAD4 corner strain output - done

        ('SDISP', 'SSDSET', 'SSDMEDIA', 'SSDFMT'),
        ('SVELO', 'SSVSET', 'SSVMEDIA', 'SSVFMT'),
        ('SACCE', 'SSASET', 'SSAMEDIA', 'SSAFMT'),

        ('GPFORCE', 'GPFSET', 'GPFMEDIA', 'GPFFMT'),
        ('ESE', 'ESESET', 'ESEMEDIA', 'ESEFMT'),
        ('AEROF', 'ARFPTSET', 'ARFMEDIA', 'ARFFMT'),
        ('SVECTOR', 'SVSET', 'SVMEDIA', 'SVFMT'),
        ('MPRES', 'FLUPTSET', 'FLUMEDIA', 'FLUFMT'),
        ('ELSDCON', 'ESDPTSET', 'ESDMEDIA', 'ESDFMT'),
        ('GPSDCON', 'GSDPTSET', 'GSDMEDIA', 'GSDFMT'),

        ('BOUTPUT', 'CNTSET', 'CNTMEDIA', 'CNTFMT'),
        ('NLSTRESS', 'NLSSET', 'NLSMEDIA', 'NLSFMT'),
        ('GPSTRAIN', 'GPEPTSET', 'GPEMEDIA', 'GPEFMT'),

        ('EKE', 'EKEPTSET', 'EKEMEDIA', 'EKEFMT'), # EKETHRSH
        ('EDE', 'EDEPTSET', 'EDEMEDIA', 'EDEFMT'), # EDETHRSH
        ('GPKE', 'GPKESET', 'GPKEMEDI', 'GPKEFMT'),#

        ('ACFPM', 'ACFPMSET', 'ACFPMMED', 'ACFPMFMT'),
        ('ERP', 'ERPSID', 'ERPMED', 'ERPFMT'),
        ('NLLOAD', 'NONPTSET', 'NONMEDIA', 'NONFMT'),
        ('GPSTRESS', 'GPSPTSET', 'GPSMEDIA', 'GPSFMT'),
        ('STRFIELD', 'STFSET', 'STFMEDIA', 'STFFMT'),

        #('RCROSS', 'RCRSET', '', 'RCRFMT'),
        ('PFGRID', 'PFGSID', 'PFGMED', 'PFGFMT'),
        ('PFMODE(STRUCTURE)', 'PFMSSID', 'PFMSMED', 'PFMSFMT'),
        ('PFMODE(FLUID)', 'PFMFSID', 'PFMFMED', 'PFMFFMT'),
        ('GVECTOR', 'GVSET', 'GVMEDIA', 'GVFMT'),
        ('INTENSITY', 'INTENSET', 'INTENMED', 'INTENFMT'),

        ('DATAREC', 'DATSET', 'DATMEDIA', 'DATFMT'),
        ('VUGRID', 'VUGSET', 'VUGMEDIA', 'VUGFMT'),
        ('MAXMIN', 'MXMNGSET', 'MXMNGMDA', 'MXMNGFMT'), # MXMNESET, MXMNEMDA, MXMNEFMT

        #('MCFRACTION', 'MCFRSET', '', ''), # MCFRSOLN, MCFRSOLN, MCFROPT
        #('FATIGUE', 'FATIGUE', 'FTGMED', ''),
        ('ICF', 'ICFSET', 'ICFMED', 'ICFFMT'), # ICFGENST, ICFGENNM, ICFUSEST, ICFUSENM
        ('ACPOWER', 'ACPOWSET', 'ACPOWMED', 'ACPOWFMT'), # ACPOWCSV
        ('MODALSE', 'MDLSSET', 'MDLSMEDIA', 'MDLSFMT'),  # MDLSTHRE, MDLSTFVL
        ('MODALKE', 'MDLKSET', 'MDLKMEDIA', 'MDLKFMT'),  # MDLKESRT, MDLKTHRE, MDLKTFVL
        #('CMSENRGY', 'CMSESET', 'CMSEMDIA', ''),  # CMSEOPTS, CMSETHRE, CMSETOPN
    ]
    string_types = [
        'SUBTITLE', 'LABEL', 'TITLE', #  trivial
        'ANALYSIS',
        #'AECONFIG', # 'POSTO2NM',
        # matrices
        'K2PP', 'M2PP', 'B2PP', 'P2G',
        'K2GG', 'M2GG', 'B2GG',
        'A2GG', 'K42GG',
    ]
    special_names = [
        'AESYMXY', 'AESYMXZ',
        'VONMISES', # vonmises flag for STRESS

        'GPQSTRS', # CQUAD4 corner stress output
        'GPQFORC', # CQUAD4 corner force output
        'GPQSTRN', # CQUAD4 corner strain output

        #'MXMNESET', # MAXMIN element set selection
        #'MXMNEMDA', # MAXMIN element output media
        #'MXMNEFMT', # MAXMIN element output format

        'MCFRSET',  # MCFRACTION: Modal contribution fraction set identification number
        'MCFRSOLN', # MCFRACTION: SOLUTION=1001
        'MCFRFILT', # MCFRACTION: FILTER=0.01
        #'MCFROPT',  # MCFRACTION: options bit pattern
    ]
    integer_mapped_names = {}
    stress_mapped_names = []
    special_names = []
    subcases = cc.subcases
    for name_str in string_types:
        assert name_str in names, name_str
        values = hdf5_case_control[name_str]
        names.remove(name_str)
        for sid, value in zip(sids, values):
            value = value.strip().decode('latin1')
            if len(value) == 0:
                continue
            options = []
            subcase: Subcase = subcases[sid]
            assert isinstance(value, str), value
            subcase.add(name_str, value, options, 'STRING-type')


    for mapped_namesi in stress_mapped_names:
        (name_str, set_str, media_str, fmt_str) = mapped_namesi
        #if name_str not in names:
            #log.error(f'skipping {name_str} because it is missing...')
            #continue
        values = hdf5_case_control[set_str]
        medias = hdf5_case_control[media_str]
        fmts = hdf5_case_control[fmt_str]
        #names.remove(name_str)
        names.remove(set_str)
        names.remove(media_str)
        names.remove(fmt_str)
        for sid, value, media, fmt in zip(sids, values, medias, fmts):
            if value == 0:
                #log.debug(f'stress-type: {name_str} is empty; value={value}')
                continue
            elif value == -1: # default
                value = 'ALL'
                #log.debug(f'stress-type: {name_str} is default; value={value}')
                #continue

            options = []
            if fmt == 0: # ??
                pass # options.append('SORT1')
            elif fmt == 1: # ??
                options.append('SORT1')
            elif fmt == 2: # ??
                options.append('SORT1')
            else:
                raise RuntimeError((name_str, sid, value, media, fmt))

            if media in (1,    3,    5,    7):
                options.append('PLOT')
            if media in (   2, 3,       6, 7):
                options.append('PRINT')
            if media in (         4, 5, 6, 7):
                options.append('PUNCH')

            #if media == 1:
                #options.append('PRINT')
            #elif media == 2:
                #options.append('PLOT')
            #elif media == 3:
                #options.extend(['PLOT', 'PRINT'])
            #elif media == 4:
                #options.append('PUNCH')
            #elif media == 5:
                #options.extend(['PRINT', 'PUNCH'])
            #elif media == 6:
                #options.extend(['PLOT', 'PUNCH'])
            #else:
                #raise RuntimeError((name_str, sid, value, media, fmt))

            subcase: Subcase = subcases[sid]
            subcase.add(name_str, value, options, 'STRESS-type')


    names_copy = deepcopy(names)
    #skip_fields = ['NONFMT', 'NONMEDIA', '']
    for name in names:
        names_copy.remove(name)
        values = hdf5_case_control[name]
        for sid, value in zip(sids, values):
            if name in integer_mapped_names:
                updated_name, default = integer_mapped_names[name]
                if value == default:
                    continue
                # SPC = 1
                subcase: Subcase = subcases[sid]
                options = []
                subcase.add(updated_name, value, options, 'STRESS-type')
            elif name in special_names:
                continue
            #elif name.endswith('FMT'):
                #log.error(f'format skipping {name}={value}')
            elif name == 'DOMAIN_ID':
                log.warning(f'skipping {name}={value}')
                break
            #else:
                #log.warning(f'skipping {name}={value}')
            #cc.dtype.names
    ###
    _load_special_names(special_names, hdf5_case_control,
                        subcases, sids, cc.log)

    #str(cc)
    geom_model.case_control_deck = cc

def _load_special_names(special_names, hdf5_case_control,
                        subcases: dict[int, Subcase],
                        sids: list[int],
                        log: SimpleLogger):
    for name_str in special_names:
        values = hdf5_case_control[name_str]
        for sid, data in zip(sids, values):
            subcase: Subcase = subcases[sid]
            data_type = 'STRESS-type'

            options = []
            if name_str in ['AESYMXY', 'AESYMXZ']:
                data_type = 'STRING-type'
                if data == 0:
                    value = 'ANTISYMMETRIC'
                    continue
                # SYMMETRIC
                # ANTISYMMETRIC
                # ASYMMETRIC
                else:
                    raise NotImplementedError((name_str, data))
                subcase.add(name_str, value, options, data_type)
                continue
            elif name_str in ['GPQSTRS', 'GPQFORC', 'GPQSTRN']:  # pre-existing
                if data == 0:  # not 100%
                    continue
                if name_str == 'GPQSTRS':
                    slot = 'STRESS'
                elif name_str == 'GPQFORC':
                    slot = 'FORCE'
                elif name_str == 'GPQSTRN':
                    slot = 'STRAIN'
                else:
                    raise RuntimeError(name_str)
                options = subcase[slot][1]
                options.append('BILIN')
                continue
            elif name_str == 'VONMISES':  # pre-existing
                if 'STRESS' not in subcase:  # not set either way; maybe it's SHEAR?
                    continue
                if data == 1:
                    value = 'VONMISES'
                elif data in {3, 7}:
                    pass ## TODO: what????
                    log.warning(f'skipping VONMISES={data}')
                    #value = 'VONMISES'
                    continue
                else:
                    raise NotImplementedError((name_str, data))
                options = subcase['STRESS'][1]
                options.append(value)
                continue
            elif name_str == 'MCFRSET':  #  new
                if data == 0:
                    continue
                subcase.add('MCFRACTION', data, options, data_type)
                continue
            elif name_str == 'MCFRSOLN':  # pre-existing
                if data == 0:
                    continue
                options = subcase['MCFRACTION'][1]
                options.append(f'SOLUTION={data}')
                continue
            elif name_str == 'MCFRFILT':  # pre-existing
                if data == 0.0:
                    continue
                options = subcase['MCFRACTION'][1]
                options.append(f'FILTER={data}')
                continue
            elif name_str == 'MCFROPT':
                # pre-existing
                #options = subcase['MCFRACTION'][1]
                # ???
                #options.append(f'FILTER={value}')
                continue
            else:
                raise NotImplementedError(name_str)
            raise RuntimeError('asdf')
    for subcase in subcases:
        str(subcase)
