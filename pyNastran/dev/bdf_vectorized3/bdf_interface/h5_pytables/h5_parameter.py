from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

from .utils import cast_encoding_strip, cast_encoding_strip0, get_group_name
from pyNastran.utils.numpy_utils import integer_types, float_types
from pyNastran.bdf.case_control_deck import CaseControlDeck
if TYPE_CHECKING: # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def _load_h5_parameter_pvt(model: BDF, group: Group, encoding: str):
    log = model.log
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            group_name = get_group_name(group)
            raise RuntimeError(group_name)
            #msg += f'{indent}Group: {group_name}\n'
            #indent += '    '
            #print_group(group, indent=indent)
            #continue
        elif class_id == 'TABLE':
            group_name = get_group_name(group)
            #log.warning(f'PARAM table name={group_name}')
        else: # pragma: no cover
            raise RuntimeError(class_id)

        if group_name == 'PVT':
            datatype = h5_element.name
            data = h5_element.read()
            name_bytes = data['NAME']
            name = cast_encoding_strip(name_bytes, encoding)
            if datatype == 'INT':
                # NAME, VALUE
                value = data['VALUE']
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            elif datatype == 'CHAR':
                # NAME, VALUE
                value_bytes = data['VALUE']
                value = cast_encoding_strip(value_bytes, encoding)
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            elif datatype in {'DOUBLE'}:  # TODO: does FLOAT exist?
                # NAME, VALUE
                value = data['VALUE']
                for key, valuei in zip(name, value):
                    model.add_param(key, [valuei])
            else: # pragma: no cover
                raise NotImplementedError(datatype)
        else: # pragma: no cover
            raise NotImplementedError(group_name)
        x = 1

def _load_h5_parameter_group(model: BDF, group: Group, encoding: str):
    log = model.log
    group_name = get_group_name(group)
    skip_names = {'TSTEP', } # , 'PVT',
    if group_name in skip_names:
        model.log.warning(f'skipping {group_name} in _load_h5_parameter_group')
    elif group_name == 'CASECC':
        subcase_data = group['SUBCASE'].read()
        #print(f'subcase_data = {subcase_data}')
        #print(f'  dtype subcase_data = {subcase_data.dtype}')
        #for h5_elementi in group._f_iter_nodes():
            #print(h5_elementi)
            #print(h5_elementi.dtype.name)
        nsubcase = len(subcase_data)
        casecc = {sid: {} for sid in subcase_data['SID']}
        cc = CaseControlDeck([])
        for sid in casecc:
            cc.create_new_subcase(sid)
        model.case_control_deck = cc
        subcases = cc.subcases

        integer_map = {
            'SPCSET': 'SPC',
            'MPCSET': 'MPC',
            'TIC': 'TIC',
            'TFSET': 'TF',
            'DESOBJ': 'DESOBJ',
            'DESSUB': 'DESSUB',
            'FLUTTER': 'FMETHOD',
        }
        matrices = {
            'M2GG', 'B2GG', 'K2PP', 'B2PP', 'K2GG', 'P2G',
            'K42GG', 'A2GG', 'M2PP', }

        sids = subcase_data['SID']
        for key in subcase_data.dtype.names:
            data = subcase_data[key]
            data2 = data.flatten()
            #log.info(f'{key}: data={data} data2={data2}')

            value = data2[0]
            if isinstance(value, integer_types) and data.max() == 0 and data.min() == 0:
                continue
            if isinstance(value, float_types) and data.max() == 0. and data.min() == 0.:
                continue
            elif isinstance(value, str):
                data = np.array([value.strip() for value in data])
                value = data[0]
                if np.all(data == ''):
                    continue
            elif isinstance(value, bytes):
                data = np.array([value.strip() for value in data])
                value = data[0]
                if np.all(data == b''):
                    continue

            if key in integer_map:
                mapped_name = integer_map[key]
                for (sid, valuei) in zip(sids, data):
                    if valuei == 0:
                        continue
                    subcases[sid].add_integer_type(mapped_name, valuei)
                log.debug(f'int {key}: {data}')
                continue
            elif key in ['SUBTITLE', 'LABEL']:  # bytes
                for (sid, valuei) in zip(sids, data):
                    svalue = valuei.decode('latin1')
                    if valuei == '':
                        continue
                    subcases[sid].add_string_type(key, svalue)
                log.debug(f'str {key}: {data}')
                continue
            elif key in matrices:  # bytes?
                for (sid, valuei) in zip(sids, data):
                    svalue = valuei.decode('latin1')
                    if valuei == '':
                        continue
                    raise TypeError(svalue)
            #log.debug(f'{key}: {data}; value={value!r} type={type(value)}')

        # has a bug...
        #log.debug(str(cc))
        #
        # ints:
        # MPCSET, SPCSET, TIC, TFSET, DESOBJ, DESSUB,
        #
        # SPCFORCE:
        # MPCFORCE: 'MPCFSET', 'MPCMEDIA', 'MPCFFMT',
        # DISP:
        # ACCEL:
        # VELO:
        # STRESS:
        # STRAIN:
        # GPKE:      'GPKESET', 'GPKEMEDI', 'GPKEFMT',
        # ESE:
        #
        # matrices:
        #  'M2GG', 'B2GG',
        # bytes:
        # SUBTITLE, LABEL,
        #
        # {'names': ['SID', 'ESLSET', 'REESET', 'ELDSET', 'THLDSET', 'THMATSET', 'NONPTSET', 'NONMEDIA', 'NONFMT', 'DYMLDSET', 'FEQRESET', , 'SYMFLG', 'LDSPTSET', 'LDSMEDIA', 'LDSFMT', 'DPLPTSET', 'DPLMEDIA', 'DPLFMT', 'STSPTSET', 'STSMEDIA', 'STSFMT', 'FCEPTSET', 'FCEMEDIA', 'FCEFMT', 'ACCPTSET', 'ACCMEDIA', 'ACCFMT', 'VELPTSET', 'VELMEDIA', 'VELFMT', 'FOCPTSET', 'FOCMEDIA', 'FOCFMT', 'TSTEPTRN', 'TITLE', 'STPLTFLG', 'AXSYMSET', 'NOHARMON', 'TSTRU', 'SETKP', 'FLAGKP', 'SETMP', 'FLAGMP', 'SETBP', 'FLAGBP', 'OUTRESPV', 'SEDR', 'FLDBNDY', 'CEESET', 'DAMPTBL', 'DYNRED', 'SSDSET', 'SSDMEDIA', 'SSDFMT', 'SSVSET', 'SSVMEDIA', 'SSVFMT', 'SSASET', 'SSAMEDIA', 'SSAFMT', 'NONLINLD', 'PARTIT', 'CYCLIC', 'RANDOM', 'NONPARAM', 'LCC', 'GPFSET', 'GPFMEDIA', 'GPFFMT', 'ESESET', 'ESEMEDIA', 'ESEFMT', 'ARFPTSET', 'ARFMEDIA', 'ARFFMT', 'SEID', 'LCN', 'GUST', 'SEFINAL', 'SEMG', 'SEKR', 'SELG', 'SELR', 'SEEX', 'SETKG', 'FLAGKG', 'SETMG', 'FLAGMG', 'SETBG', 'FLAGBG', 'SVSET', 'SVMEDIA', 'SVFMT', 'FLUPTSET', 'FLUMEDIA', 'FLUFMT', 'HOUT', 'NOUT', 'SETPG', 'FLAGPG', 'LOADSET', 'SEMR', 'VONMISES', 'SECMDFLG', 'GPSPTSET', 'GPSMEDIA', 'GPSFMT', 'STFSET', 'STFMEDIA', 'STFFMT', 'CLOAD', 'SET2ID', 'DSAPRT', 'DSASTORE', 'DSAOUTPT', 'STNSET', 'STNMEDIA', 'STNFMT', 'APRESS', 'TRIM', 'MODLIST', 'REESETF', 'ESDPTSET', 'ESDMEDIA', 'ESDFMT', 'GSDPTSET', 'GSDMEDIA', 'GSDFMT', 'SEDV', 'SERE', 'SERS', 'CNTSET', 'CNTMEDIA', 'CNTFMT', 'DIVERG', 'OUTRCV', 'STATSUBP', 'DFT1', 'DFT2', 'ADAPT', 'SUBSPAN', 'DESGLB', 'ANALYSIS', 'GPQSTRS', 'GPQFORC', 'GPQSTRN', 'SUPORT1', 'STATSUBB', 'BCID', 'AUXMODEL', 'ADACT', 'DATSET', 'DATMEDIA', 'DATFMT', 'VUGSET', 'VUGMEDIA', 'VUGFMT', 'REUESET', 'DAMPTBLF', 'ITERMETH', 'NLSSET', 'NLSMEDIA', 'NLSFMT', 'MODTRKID', 'DSAFORM', 'DSAEXPO', 'DSABEGIN', 'DSAINTVL', 'DSAFINAL', 'DSASETID', 'SORTFLG', 'RANDBIT', 'AECONFIG', 'AESYMXY', 'AESYMXZ', 'OCIDREQ', 'GPEPTSET', 'GPEMEDIA', 'GPEFMT', 'TEMPMAT', 'AECSSSET', 'EKEPTSET', 'EKEMEDIA', 'EKEFMT', 'EKETHRSH', 'EDEPTSET', 'EDEMEDIA', 'EDEFMT', 'EDETHRSH', 'DFT3', 'DFR1', 'DFR2', 'DFR3', 'SETK4G', 'FLAGK4G', 'NK42GG', 'NA2GG', 'EFFMASET', 'EFFMAGID', 'EFFMATHR', 'EQUILMED', 'EQUILGRD', 'RCRSET', 'RCRFMT', 'AEUXREF', 'GCHK', 'GCHKOUT', 'GCHKSET', 'GCHKGID', 'GCHKTHR', 'GCHKRTHR', 'GCHKDREC', 'ASPCMED', 'ASPCEPS', 'ASPCPRT', 'ASPCPCH', 'NK2PP', 'NM2PP', 'NB2PP', 'NK2GG', 'NM2GG', 'NB2GG', 'NP2G', 'GEODSET', 'GEODMXMN', 'GEODOCID', 'GEODNUMB', 'GEOLSET', 'GEOLMXMN', 'GEOLOCID', 'GEOLNUMB', 'GEOSSET', 'GEOSMXMN', 'GEOSOCID', 'GEOSNUMB', 'GEOMSET', 'GEOMMXMN', 'GEOMOCID', 'GEOMNUMB', 'GEOASET', 'GEOAMXMN', 'GEOAOCID', 'GEOANUMB', 'GEOVSET', 'GEOVMXMN', 'GEOVOCID', 'GEOVNUMB', 'NTFL', 'SEDAMP', 'WCHK', 'WCHKOUT', 'WCHKSET', 'WCHKGID', 'WCHKCGI', 'WCHKWM', 'EXSEOUT', 'EXSEMED', 'EXSEUNIT', 'EXSERES1', 'EXSERES2', 'FK2PP', 'FM2PP', 'FB2PP', 'FK2GG', 'FM2GG', 'FB2GG', 'TICTYPE', 'FK42GG', 'FA2GG', 'SUBSTEP', 'STEPID', 'NSMID', 'ROUTDISP', 'ROUTVELO', 'ROUTACCE', 'ROUTLOAD', 'ROUTSPCF', 'ROUTSTRS', 'ROUTFORC', 'ROUTSTRN', 'ROUTMSCF', 'MDLSSET', 'MDLSMEDIA', 'MDLSFMT', 'MDLSESRT', 'MDLSTHRE', 'MDLSTFVL', 'MDLKSET', 'MDLKMEDIA', 'MDLKFMT', 'MDLKESRT', 'MDLKTHRE', 'MDLKTFVL', 'ACPOWSET', 'ACPOWMED', 'ACPOWFMT', 'ACPOWCSV', 'NLOUT', 'SEEFMNO', 'SEEFMHV', 'SEEFMDLF', 'SEEFMBND', 'SEEFMSDP', 'CONNECTOR', 'ESETHRSH', 'POSTUNIT', 'POSTOPT1', 'POSTOPT2', 'TICDIFF', 'DSAESEID', 'MXMNGSET', 'MXMNGMDA', 'MXMNGFMT', 'MXMNESET', 'MXMNEMDA', 'MXMNEFMT', 'MCFRSET', 'MCFRSOLN', 'MCFRFILT', 'MCFROPT', 'ELSUMID', 'ELSUMOPT', 'ELSUMDUM', 'RGYRO', 'CMSESET', 'CMSEMDIA', 'CMSEOPTS', 'CMSETHRE', 'CMSETOPN', 'GPRSORT', 'MASSSET', 'AESOLN', 'POSTO2NM', 'RANDVAR', 'RSVCRQTS', 'RSVCOPTS', 'RSVCSTBS', 'RSVCRQTC', 'RSVCOPTC', 'RSVCSTBC', 'DESVAR', 'BCONTACTI', 'BCONTACTC', 'MODSELS1', 'MODSELS2', 'MODSELS3', 'MODSELS4', 'MODSELS5', 'MODSELS6', 'MODSELS7', 'MODSELS8', 'MODSELS9', 'MODSELF1', 'MODSELF2', 'MODSELF3', 'MODSELF4', 'MODSELF5', 'MODSELF6', 'MODSELF7', 'MODSELF8', 'MODSELF9', 'FTNURN', 'SUFNAM1', 'SUFNAM2', 'ENVELOP1', 'ENVELOP2', 'GPFLXSET', 'GPFLXMED', 'CAMPBELL', 'SPLINOUT', 'MONITOR', 'FBODYLD', 'STOCHAST', 'EXPTLDID', 'EXPTLDNM', 'EXPTLDSI', 'AERCONFIG', 'NLHARM', 'PFMSSID', 'PFMSMED', 'PFMSFMT', 'PFMSSSET', 'PFMSFLTR', 'PFMSOPTS', 'PFMSSTMP', 'PFMFSID', 'PFMFMED', 'PFMFFMT', 'PFMFSSET', 'PFMFPSET', 'PFMFFLTR', 'PFMFOPTS', 'PFMFFLMP', 'PFMFSTMP', 'PFPSID', 'PFPMED', 'PFPFMT', 'PFPSSET', 'PFPFLTR', 'PFPOPTS', 'PFPPSET', 'PFGSID', 'PFGMED', 'PFGFMT', 'PFGSSET', 'PFGGSET', 'NLICCASE', 'NLICSTEP', 'NLICLFAC', 'ACFPMSET', 'ACFPMMED', 'ACFPMFMT', 'FRFFLAG', 'FRFCMPID', 'FRFCONST', 'FRFUNTNO', 'FRFCMPNM', 'DESMOD', 'RSDAMPST', 'RSDAMPFL', 'VCCT', 'FP2G', 'FRQVAR', 'HADAPT', 'BCHANGE', 'BCMOVE', 'BSQUEAL', 'UNGLUE', 'HSUBCASE', 'HSTEP', 'HTIME', 'ERPSID', 'ERPMED', 'ERPFMT', 'ERPSSET', 'ERPFLTR', 'ERPOPTS', 'ERPCSV', 'TESTTHRR', 'TESTTHRI', 'ASMOUTFL', 'ASMOUTNM', 'ICFUNTNO', 'TFOSET', 'TFUNIT', 'TFLSET', 'TFOPTS', 'NLOOPH', 'NLSTEP', 'NMODES', 'RCPARM', 'NLOPCTRL', 'NLOPDBG', 'NLOPPOST', 'NLOPMPCH', 'NLICSSTP', 'DEACTEL', 'ACTIVAT', 'INTENSET', 'INTENMED', 'INTENFMT', 'NLICTOLR', 'PFPSSID', 'PFPSMED', 'PFPSFMT', 'PFPSSSET', 'PFPSFLTR', 'PFPSOPTS', 'PFPSPSET', 'ACTISET', 'ELSOSET', 'ELSRSET', 'ELSSSET', 'ELSTHRS', 'ELSBITS', 'WTSOSET', 'WTSRSET', 'WTSSSET', 'WTSTHRS', 'WTSBITS', 'PACOSET', 'PACSSET', 'PACBITS', 'IRLOAD', 'ICFSET', 'ICFMED', 'ICFFMT', 'ICFGENST', 'ICFGENNM', 'ICFUSEST', 'ICFUSENM', 'HISTSET', 'HISTTYPE', 'HISTFMT', 'FATIGUE', 'FTGMED', 'FTGFMT', 'NLOPDELI', 'NLOPGRID', 'DASAVE', 'GVSET', 'GVMEDIA', 'GVFMT', 'ERMPF', 'NVELOSET', 'NVELFDEF', 'NVELFDTA', 'NVELTHRS', 'NVELBITS', 'VITSOSET', 'VITSRSET', 'VITSSSET', 'VITSTHRS', 'VITSBITS', 'THLDVER', 'THMATVER', 'NPEAK', 'NEAR', 'LFREQ', 'HFREQ', 'RTYPE', 'PSCALE', 'NSAMP', 'MONSET', 'SEED', 'OFFD', 'FSORT2', 'MFREQ', 'BCONTACT', 'SOLNID', 'ExtDROut', 'ExtDROmd', 'ExtDROun', 'ExtDRIn', 'ExtDRImd', 'ExtDRIun', 'FemCheck', 'LSEM', 'SYM_LEN', 'SYM_POS', 'DOMAIN_ID']
        #datatype = group.name
        #model.log.warning(f'skipping {group_name}={datatype} in _load_h5_parameter_group')
        #data = group.read()
        #log.warning(type(group))
    elif group_name == 'PVT':
        _load_h5_parameter_pvt(model, group, encoding)
    elif group_name == 'AECOMP':
        #('NAME', 'LISTTYPE', 'LISTIDS_POS', 'LISTIDS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        list_type = cast_encoding_strip(identity['LISTTYPE'], encoding)
        nlists = identity['LISTIDS_LEN']
        lists = group['LISTIDS'].read()
        model.aecomp._save(name, list_type, nlists, lists)
    elif group_name == 'AECOMPL':
        #('NAME', 'LABELS_POS', 'LABELS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        name = identity['NAME']
        nlabels = identity['LABELS_LEN']
        labels_bytes = group['LABELS'].read()
        labels = cast_encoding_strip0(labels_bytes, encoding)
        model.aecompl._save(name, nlabels, labels)
    elif group_name == 'AEFACT':
        #('SID', 'DIS_POS', 'DIS_LEN', 'DOMAIN_ID')
        identity = group['IDENTITY'].read()
        aefact_id = identity['SID']
        nfractions = identity['DIS_LEN']
        fractions = group['DIS'].read()
        model.aefact._save(aefact_id, nfractions, fractions)
    else:
        raise NotImplementedError(group_name)

def load_h5_parameter(model: BDF, group: Group, encoding: str='latin1'):
    for h5_element in group._f_iter_nodes():
        class_id = h5_element._c_classid
        if class_id == 'GROUP':
            _load_h5_parameter_group(model, h5_element, encoding)
            continue
        elif class_id == 'TABLE':
            pass
        else: # pragma: no cover
            raise RuntimeError(class_id)

        name = h5_element.name
        data = h5_element.read()
        #if name in {'PVT'}:
            #model.log.warning(f'skipping {element_name} in load_h5_parameter')
        if name == 'TSTEPNL':
            model.log.warning(f'skipping {name} in load_h5_parameter')
        elif name == 'AEROS':
            #('ACSID', 'RCSID', 'REFC', 'REFB', 'REFS', 'SYMXZ', 'SYMXY', 'DOMAIN_ID')
            acsid = data['ACSID'][0]
            rcsid = data['RCSID'][0]
            bref = data['REFB'][0]
            cref = data['REFC'][0]
            sref = data['REFS'][0]
            sym_xz = data['SYMXZ'][0]
            sym_xy = data['SYMXY'][0]
            aeros = model.add_aeros(cref, bref, sref, acsid=acsid,
                                    rcsid=rcsid, sym_xz=sym_xz, sym_xy=sym_xy, comment='')
            str(aeros)
        elif name == 'AESTAT':
            #('ID', 'LABEL', 'DOMAIN_ID')
            aestat_ids = data['ID']
            labels = data['LABEL']
            for aestat_id, label in zip(aestat_ids, labels):
                label_str = label.strip().decode('latin1')
                model.add_aestat(aestat_id, label_str, comment='')
        elif name == 'MDLPRM':
            #('NAME', 'VALUE', 'DOMAIN_ID')
            names = data['NAME']
            values = data['VALUE']
            mdlprm_dict = {}
            for name, value in zip(names, values):
                name_str = name.strip().decode('latin1')
                mdlprm_dict[name_str] = value
            assert len(mdlprm_dict) == 1, mdlprm_dict
            assert 'HDF5' in mdlprm_dict, mdlprm_dict
            model.add_mdlprm(mdlprm_dict, comment='')
        elif name == 'MONPNT1':
            name = data['NAME']
            label = data['LABEL']
            axes = data['AXES']
            comp = data['COMP']
            cp = data['CP']
            xyz = np.stack(
                [data['X'], data['Y'], data['Z']],
                axis=1)
            assert xyz.shape[1] == 3, xyz.shape
            cd = data['CD']
            model.monpnt1._save(
                name, label, comp, axes, cp, xyz, cd)
        else: # pragma: no cover
            raise NotImplementedError(name)
