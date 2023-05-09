from typing import Optional, Any, cast
import numpy as np
from cpylog import SimpleLogger
from pyNastran.converters.tecplot.zone import CaseInsensitiveDict, Zone


def read_header_lines(lines: list[str], iline: int, line: str,
                      log: SimpleLogger,
                      # iline, title_line, header_lines, line
                      ) -> tuple[int, str, list[str], str]:
    """
    reads a tecplot header

    Examples
    --------
    **Example 1**

    TITLE     = "tecplot geometry and solution file"
    VARIABLES = "x"
    "y"
    "z"
    "rho"
    "u"
    "v"
    "w"
    "p"
    ZONE T="\"processor 1\""
    n=522437, e=1000503, ZONETYPE=FEBrick
    DATAPACKING=BLOCK

    **Example 2**

    title="Force and Momment Data for forces"
    variables="Iteration"
    "C_L","C_D","C_M_x","C_M_y","C_M_z""C_x","C_y","C_z","C_Lp","C_Dp"
    "C_Lv", "C_Dv","C_M_xp", "C_M_yp","C_M_zp"
    "C_M_xv","C_M_yv""C_M_zv","C_xp","C_yp","C_zp","C_xv","C_yv""C_zv
    "Mass flow","<greek>r</greek>","u"
    "p/p<sub>0</sub>","T","p<sub>t</sub>/p<sub>0</sub>"
    "T<sub>t</sub>","Mach"
    "Simulation Time"
    zone,t="forces"
    """
    i = 0
    title_line = ''
    #variables_line = ''
    active_key = None

    vars_found = []
    header_lines = []
    #print('-----------------------------')
    #for iii, linei in enumerate(lines):
        #if iii > 10:
            #break
        #print(linei)
    #print('-----------------------------')
    while i < 30:
        #print(iline, i, line.strip())
        #self.n = 0
        #log.info(f'{iline:d}: {line}')
        if len(line) == 0 or line[0] == '#':
            line = lines[iline].strip()
            iline += 1
            i += 1
            continue
        if line[0].isdigit() or line[0] == '-':
            #print(line)
            log.debug('breaking after finding header lines...')
            break

        uline = line.upper()
        uline2 = uline.replace(' ', '')
        if 'TITLE=' in uline2:
            title_line += line
            vars_found.append('TITLE')
            active_key = 'TITLE'
        elif 'VARIABLES' in uline2:
            vars_found.append('VARIABLES')
            #variables_line += line
            active_key = 'VARIABLES'
        else:
            #if 'ZONE T' in line:
                #vars_found.append('ZONE T')
            if 'ZONE' in uline2:
                vars_found.append('ZONE')
                active_key = 'ZONE'
            #if 'ZONE N' in uline:
                #vars_found.append('N')
            if 'ZONETYPE' in uline2:
                vars_found.append('ZONETYPE')
                active_key = 'ZONE'
            if 'DATAPACKING' in uline2:
                vars_found.append('DATAPACKING')
                active_key = 'ZONE'

        #print(active_key, line)
        if active_key in ['ZONE', 'VARIABLES']:
            header_lines.append(line.strip())
        #if len(vars_found) == 5:
            #break

        #if active_key
        i += 1
        line = lines[iline].strip()
        #log.info(f'*{iline:d}: {line}')
        iline += 1

    log.debug('vars_found = %s' % vars_found)
    #print('header_lines', header_lines)
    #print("title = %r" % title_line)
    #print("variables_line = %r" % variables_line)

    return iline, title_line, header_lines, line

def header_lines_to_header_dict(title_line: str, header_lines: list[str],
                                variables: list[str], log: SimpleLogger):
    """parses the parsed header lines"""
    #print('header_lines', header_lines)
    #headers_dict = {}
    headers_dict = CaseInsensitiveDict()
    if title_line:
        title_sline = title_line.split('=', 1)
        title = title_sline[1]
    else:
        title = 'tecplot geometry and solution file'
    headers_dict['TITLE'] = title

    if len(header_lines) == 0:
        #raise RuntimeError(header_lines)
        return None

    header = _join_headers(header_lines)

    # this is so overly complicataed and probably not even enough...
    # what about the following 'quote' style?
    headers = split_headers(header, log)
    #headers = header.replace('""', '","').split(',')

    #TITLE = "Weights=1/6,6,1"
    #Variables = "x","y","z","psi"
    #Zone N = 125, E = 64, DATAPACKING = POINT, ZONETYPE = FEBRICK

    nheaders = len(headers) - 1
    for iheader, header in enumerate(headers):
        header = header.strip()
        #print(f'{iheader} {header!r}')

    for iheader, header in enumerate(headers):
        header = header.strip()
        #print('%2i %s' % (iheader, header))
        #print('iheader=%s header=%r' % (iheader, header))
        if '=' in header:
            sline = header.split('=', 1)
            parse = False
            #print('iheader=%s nheaders=%s' % (iheader, nheaders))
            if iheader == nheaders:
                parse = True
            elif '=' in headers[iheader + 1]:
                parse = True
        elif header.upper() == 'ZONE':
            # apparently the null key is also a thing...
            # we'll use 'ZONE' because...
            headers_dict['ZONE'] = None
            parse = True
            #continue
        elif '"' in header:
            sline += [header]
            parse = False
            if iheader == nheaders:
                parse = True
            elif '=' in headers[iheader + 1]:
                parse = True
        else:
            raise NotImplementedError('header=%r headers=%r' % (header, headers))

        if parse:
            # ZONE T="FUSELAGE" I=21 J=49 K=1 F=BLOCK
            #print('  parsing')
            log.debug(f'sline = {sline}')
            key = sline[0].strip()
            ukey = key.upper()
            if ukey.startswith('ZONE '):
                # the key is not "ZONE T" or "ZONE E"
                # ZONE is a flag, T is title, E is number of elements
                ukey = ukey[5:].strip()
                headers_dict['NAME'] = ''.join(sline[1:]).strip('\"\\')

            value = [val.strip() for val in sline[1:]]
            if len(value) == 1:
                value = value[0].strip()
            #assert not isinstance(value, list), value

            if ukey == 'NODES':
                ukey = 'N'
            if ukey == 'ELEMENTS':
                ukey = 'E'
            if ukey == 'F':
                ukey = 'DATAPACKING'

            if 'DATASETAUXDATA' in ukey:
                #DATASETAUXDATA Common.AngleOfAttack="0.000000"
                #DATASETAUXDATA Common.DensityVar="5"
                #DATASETAUXDATA Common.GasConstant="0.7142857143"
                #DATASETAUXDATA Common.PressureVar="9"
                #DATASETAUXDATA Common.ReferenceMachNumber="1.400000"
                #DATASETAUXDATA Common.UVar="6"
                #DATASETAUXDATA Common.VectorVarsAreVelocity="TRUE"
                #DATASETAUXDATA Common.VVar="7"
                #DATASETAUXDATA Common.WVar="8"
                base, new_key = key.split('.', 1)
                assert base == 'DATASETAUXDATA Common', f'base={base!r}'
                allowed_keys = ['AngleOfAttack', 'GasConstant', 'ReferenceMachNumber', 'VectorVarsAreVelocity',
                                'DensityVar', 'PressureVar', 'UVar', 'VVar', 'WVar', ]
                param = headers_dict.get('COMMON_PARAMS', {})
                param[new_key] = value
                headers_dict['COMMON_PARAMS'] = param
                assert new_key in allowed_keys, 'new_key=%r; allowed=[%s]' % (new_key, ', '.join(allowed_keys))
            else:
                headers_dict[ukey] = value
                #print('  ', value)
                #value = value.strip()

                # 'T', 'ZONE T',  ???
                #   'DT', 'SOLUTIONTIME', 'STRANDID', # tecplot 360 specific things not supported
                allowed_keys = ['VARIABLES', 'T', 'ZONETYPE', 'DATAPACKING', # 'TITLE',
                                'N', 'E', 'F', 'DT', 'SOLUTIONTIME', 'STRANDID',
                                'I', 'J', 'K',]
                assert ukey in allowed_keys, 'ukey=%r; allowed=[%s]' % (ukey, ', '.join(allowed_keys))
            parse = False
    log.debug(f'headers_dict = {headers_dict}')
    #print(headers_dict.keys())
    #print(headers_dict)

    _simplify_header(headers_dict, variables)
    assert len(headers_dict) > 0, headers_dict
    return headers_dict

def split_headers(headers_in: str, log: SimpleLogger) -> list[str]:
    log.debug(f'headers_in = {headers_in}')
    #allowed_keys = ['TITLE', 'VARIABLES', 'T', 'ZONETYPE', 'DATAPACKING',
                    #'N', 'E', 'F', 'DT', 'SOLUTIONTIME', 'STRANDID',
                    #'I', 'J', 'K'
                    #]
    #print(f'header1 = {headers_in}')
    header = headers_in.replace('""', '","')
    #print(f'header2 = {header}')
    cheaders = header.split(',')

    #print(header)
    #print(cheaders)
    #header = cheaders[0]
    #headers = [header]
    #i = 1
    #while i < len(cheaders):
        #headeri = cheaders[i]
        #uheaderi = headeri.upper().replace(' ', '')
        #is_key = [uheaderi.startswith(key+'=') for key in allowed_keys]
        #if any(is_key):
            #print('key!', headeri)
            #header = headeri
            #headers.append(header.lstrip())
        #else:
            #headers[-1] += ',' + headeri
        #i += 1
    #print('headers')
    #for headeri in headers:
        #print('  ', headeri)

    #print(headers)
    #print(header.replace('""', '","'))
    #if ''
    #headers = header.replace('""', '","').split(',')
    return cheaders

def _join_headers(header_lines: list[str]) -> str:
    """smart join by commas"""
    header = ','.join([headeri.strip(', ') for headeri in header_lines])
    return header

def _simplify_header(headers_dict: dict[str, Any], variables: list[str]) -> None:
    """cast the integer headers and sets the variables"""
    # unstructured
    if 'N' in headers_dict: # nnodes
        headers_dict['N'] = int(headers_dict['N'])
    if 'E' in headers_dict: # nelements
        headers_dict['E'] = int(headers_dict['E'])

    # structured
    if 'I' in headers_dict:
        headers_dict['I'] = int(headers_dict['I'])
    if 'J' in headers_dict:
        headers_dict['J'] = int(headers_dict['J'])
    if 'K' in headers_dict:
        headers_dict['K'] = int(headers_dict['K'])

    #print('_simplify_header', variables, headers_dict)
    if 'TITLE' not  in headers_dict:
        headers_dict['TITLE'] = 'tecplot geometry and solution file'

    if 'VARIABLES' in headers_dict and variables is None:
        #print('VARIABLES' in headers_dict, variables is None)
        _simplify_variables(headers_dict)
    elif 'VARIABLES' in headers_dict:
        _simplify_variables(headers_dict)
    elif variables is not None:
        headers_dict['VARIABLES'] = variables
    else:
        raise RuntimeError('no variables...')

def _simplify_variables(headers_dict) -> None:
    variables = headers_dict['VARIABLES']
    headers_dict['VARIABLES'] = [var.strip('"') for var in variables]

def read_zonetype(log: SimpleLogger,
                  zone: Zone, zone_type: str,
                  lines: list[str], iline: int,
                  iblock: int,
                  headers_dict: dict[str, Any],
                  line: str,
                  nnodes: int, nelements: int,
                  zone_data_list: list[np.ndarray],
                  hexas_list: list[np.ndarray],
                  tets_list: list[np.ndarray],
                  quads_list: list[np.ndarray],
                  tris_list: list[np.ndarray],
                  data_packing: Optional[str]=None,
                  fe: Optional[str]=None) -> int:
    """
    Parameters
    ----------
    zone_type : str

    fe : str
      - a zone_type.upper() string???
      - FEPOINT

    reads:
      - ZONE E
      - ZONE T

    ZONE is a flag, T is title, E is number of elements

    -------------  ---------  ----------  ----------------------------------------------
    Parameter      Ordered    Finite      Description
                   Data       Element
    -------------  ---------  ----------  ----------------------------------------------
    T="title"      Yes        Yes         Zone title.
    I=imax         Yes        No          Number of points in 1st dimension.
    J=jmax         Yes        No          Number of points in 2nd dimension.
    K=kmax         Yes        No          Number of points in 3rd dimension.
    C=colour       Yes        Yes         Colour from WHITE, BLACK, RED, GREEN,
                                                      BLUE, CYAN, YELLOW, PURPLE,
                                                      CUST1, CUST2,....CUST8.
    F=format       Yes        Yes         POINT or BLOCK for ordered data.
                                          FEPOINT or FEBLOCK for finite element.
    D=(list)       Yes        Yes         A list of variable names to to include
                                          from the last zone.
    DT=(list)      Yes        Yes         A list of datatypes for each variable.
                                          SINGLE, DOUBLE, LONGINT, SHORTINT, BYTE, BIT.
    N=num          No         Yes         Number of nodes.
    E=num          No         Yes         Number of elements.
    ET=type        No         Yes         Element type from TRIANGLE, BRICK,
                                                            QUADRILATERAL, TETRAHEDRON.
    NV=variable    No         Yes         Variable for node value.
    -------------  ---------  ----------  ----------------------------------------------
    http://paulbourke.net/dataformats/tp/
    """
    #print('self.variables', self.variables)
    #ndim = zone.ndim
    #print('iblock =', iblock)
    if iblock == 0:
        variables = headers_dict['VARIABLES']
        variables = [variable.strip(' \r\n\t"\'') for variable in variables]
        zone.variables = variables
        log.debug('zone.variables = %s' % zone.variables)
        nresults = len(variables) # x, y, z, rho, u, v, w, p
        log.debug('nresults = %s' % nresults)

    log.debug(str(headers_dict))
    is_unstructured = False
    is_structured = False
    if zone_type in {'FETRIANGLE', 'FEQUADRILATERAL', 'FETETRAHEDRON', 'FEBRICK'}:
        #print('headers_dict =', headers_dict)
        nnodesi = headers_dict['N']
        nelementsi = headers_dict['E']
        is_unstructured = True
    elif zone_type in {'POINT', 'BLOCK'}: #  structured
        ni = headers_dict['I']
        if 'J' in headers_dict:
            nj = headers_dict['J']
            if 'K' in headers_dict:
                # 3d
                nk = headers_dict['K']
                nnodesi = ni * nj * nk
                nelementsi = (ni - 1) * (nj - 1) * (nk - 1)
            else:
                # 2d
                nnodesi = ni * nj
                nelementsi = (ni - 1) * (nj - 1)
        else:
            assert 'K' not in headers_dict, list(headers_dict.keys())
            nnodesi = ni
            nelementsi = (ni - 1)
        assert nelementsi >= 0, nelementsi
        #nelementsi = 0
        elements = None # np.zeros((nelementsi, 8), dtype='int32')
        is_structured = True
    else:
        raise NotImplementedError('zone_type = %r' % zone_type)
    log.info(f'zone_type={zone_type} data_packing={data_packing} '
             f'nnodes={nnodesi} nelements={nelementsi}')

    assert nnodesi > 0, nnodesi
    assert nresults >= 0, 'nresults=%s' % nresults

    xyz_results = np.zeros((nnodesi, nresults), dtype='float32')
    if zone_type == 'FEBRICK':
        # hex
        elements = np.zeros((nelementsi, 8), dtype='int32')
    elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETETRAHEDRON'):
        # quads / tets
        elements = np.zeros((nelementsi, 4), dtype='int32')
    elif zone_type == 'FETRIANGLE':
        # tris
        elements = np.zeros((nelementsi, 3), dtype='int32')
    #elif zone_type == 'FEBLOCK':
        #pass
    elif  zone_type in ['POINT', 'BLOCK']:
        # already handled
        #print('data')
        pass
    else:
        #if isinstance(zone_type, list):
            #raise NotImplementedError(zone_type[0])
        raise NotImplementedError(zone_type)

    sline = split_line(line.strip())
    if zone_type in ('FEBRICK', 'FETETRAHEDRON'):
        if data_packing == 'POINT':
            for inode in range(nnodesi):
                if inode == 0:
                    log.debug('zone_type=%s sline=%s' %(zone_type, sline))
                #if not len(sline[3:]) == len(results[inode, :]):
                    #msg = 'sline[3:]=%s results[inode, :]=%s' % (sline[:3], results[inode, :])
                    #raise RuntimeError(msg)

                try:
                    xyz_results[inode, :] = sline
                except ValueError:
                    msg = 'i=%s line=%r\n' % (inode, line)
                    msg += 'sline = %s' % str(sline)
                    print(msg)
                    raise
                iline, line, sline = get_next_sline(lines, iline)
        elif data_packing == 'BLOCK':
            iline, line, sline = read_zone_block(
                lines, iline, xyz_results, nresults, zone_type,
                sline, nnodesi, log)
            #print('sline =', sline)
        else:
            raise NotImplementedError(data_packing)
    elif zone_type in ('FEPOINT', 'FEQUADRILATERAL', 'FETRIANGLE'):
        iline, line, sline = _read_zonetype_fe(
            iline, line, lines, nnodesi, xyz_results)

    elif zone_type == 'POINT':
        nvars = len(zone.variables)
        iline, line, sline = read_point(
            lines, iline, xyz_results, zone_type,
            line, sline, nnodesi, nvars, log)
    elif zone_type == 'BLOCK':
        nvars = len(zone.variables)
        iline, line, sline = read_block(
            lines, iline, xyz_results, zone_type,
            line, sline, nnodesi, nvars, log)
    else:  # pragma: no cover
        raise NotImplementedError(zone_type)

    #print(elements.shape)
    #print('xyz[0 , :]', xyz[0, :])
    #print('xyz[-1, :]', xyz[-1, :])
    #print(sline)
    if is_structured:
        pass
    elif is_unstructured:
        iline, line, sline = read_unstructured_elements(
            lines, iline, sline, elements, nelementsi)

        #print(f.readline())

        if zone_type == 'FEBRICK':
            hexas_list.append(elements + nnodes)
        elif zone_type == 'FETETRAHEDRON':
            tets_list.append(elements + nnodes)
        elif zone_type in ('FEPOINT', 'FEQUADRILATERAL'):
            # TODO: why are points stuck in the quads?
            quads_list.append(elements + nnodes)
        elif zone_type == 'FETRIANGLE':
            tris_list.append(elements + nnodes)
        else:
            raise NotImplementedError(zone_type)
    else:
        raise RuntimeError()
    zone_data_list.append(xyz_results)
    nnodes += nnodesi
    nelements += nelementsi
    log.debug('nnodes=%s nelements=%s (0-based)' % (nnodes, nelements))
    del headers_dict
    iblock += 1
    if iblock == 10:
        return
    log.debug('final sline=%s' % sline)
    return iline

def _read_zonetype_fe(iline: int,
                      line: str,
                      lines: list[str],
                      nnodesi: int,
                      xyz_result: np.ndarray) -> tuple[int, str, list[str]]:
    """
    reads:
     - FEPOINT
     - FEQUADRILATERAL
     - FETRIANGLE
    """
    sline = split_line(line.strip())
    nexpected = xyz_result.shape[1]
    if len(sline) == nexpected:
        for inode in range(nnodesi):
            #print(iline, inode, sline, len(sline))
            xyz_result[inode, :] = sline
            #if abs(xyz[inode, 1]) <= 5.0:
                #msg = 'inode=%s xyz=%s'  % (inode, xyz[inode, :])
                #raise RuntimeError(msg)

            iline, line, sline = get_next_sline(lines, iline)
    else:
        nvalues_to_read = nexpected * nnodesi
        inode = 0
        #print(sline)
        while len(sline) <= nvalues_to_read:
            iline, line, slinei = get_next_sline(lines, iline)
            #print(iline, inode, slinei, len(sline))
            #if iline % 10 == 0:
                #print(iline, len(sline), nvalues_to_read)
            sline.extend(slinei)
            inode += 1
        # remove the last sline
        nslinei = len(slinei)
        xyz_result_list = sline[:-nslinei]
        assert len(xyz_result) == nvalues_to_read, (len(xyz_result_list), nvalues_to_read)
        xyz_result_temp = np.array(xyz_result_list, dtype='float32')
        xyz_result[:, :] = xyz_result_temp.reshape((nexpected, nnodesi)).T
        del sline, inode

        #print(iline, slinei)
        sline = slinei

    return iline, line, sline

def read_zone_block(lines: list[str], iline: int,
                    xyz_results: np.ndarray, nresults: int,
                    zone_type: str,
                    sline: list[str], nnodes: int,
                    log: SimpleLogger,
                    ) -> tuple[int, str, list[str]]: # iline, line, sline
    """a zone can be structured or unstructred"""
    #print('***', iline, sline)

    # read all data
    #result = sline
    #iresult = len(sline)
    #nresult = len(sline)

    xyz_result_list: list[str] = []
    iresult = 0
    nresult = 0
    nnodes_max = nresults * nnodes
    #print('nnodes_max =', nnodes_max)
    while nresult < nnodes_max: #  changed from iresult to nresult
        #print('zb', iline, sline, len(sline))
        xyz_result_list += sline
        nresult += len(sline)
        if iresult >= nnodes_max:
            log.debug('breaking...')
            #break
        iline, line, sline = get_next_sline(lines, iline)
        if iresult == 0:
            log.debug('zone_type=%s sline=%s' % (zone_type, sline))
        iresult += len(sline)
        #print('len', iresult, nresult, len(result))
    #print(result, len(result))
    #for i, value in enumerate(xyz_result_list):
        #assert '.' in value, 'i=%i value=%s' % (i, value)
    assert len(xyz_result_list) == nnodes_max, 'len(xyz_result_list)=%s expected=%s' % (len(xyz_result_list), nnodes_max)
    #-----------------

    # pack data
    for ires in range(nresults):
        i0 = ires * nnodes
        i1 = (ires + 1) * nnodes #+ 1
        if len(xyz_result_list[i0:i1]) != nnodes:
            msg = 'ires=%s len=%s nnodes=%s' % (
                ires, len(xyz_result_list[i0:i1]), nnodes)
            raise RuntimeError(msg)
        xyz_results[:, ires] = xyz_result_list[i0:i1]

    # setup
    #iline, line, sline = get_next_sline(lines, iline)
    return iline, line, sline

def read_unstructured_elements(lines: list[str], iline: int, sline: list[str],
                               elements: np.ndarray, nelements: int,
                               ) -> tuple[int, str, list[str]]: # iline, line, sline
    assert '.' not in sline[0], sline

    i = 0
    #print('nelements =', nelements)
    for i in range(nelements):
        #print(iline, i, sline)
        try:
            elements[i, :] = sline
        except IndexError:
            raise RuntimeError('i=%s sline=%s' % (i, str(sline)))
        except ValueError:
            raise RuntimeError('i=%s sline=%s' % (i, str(sline)))

        iline, line, sline = get_next_sline(lines, iline)
        #line = lines.readline()
        #iline += 1
        #sline = line.strip().split()
    return iline, line, sline


def read_point(lines: list[str], iline: int,
               xyz_results: np.ndarray,
               zone_type: str,
               line: str, sline: list[str],
               nnodes: int, nvars: int,
               log: SimpleLogger) -> tuple[int, str, list[str]]: # (iline, line, sline)
    """a POINT grid is a structured grid"""
    log.debug(f'start of POINT (structured); nnodes={nnodes} '
              f'nvars={nvars} zone_type={zone_type}')
    for inode in range(nnodes):
        iline, sline = get_next_nsline(lines, iline, sline, nvars)
        #print(iline, inode, sline)

        #if inode == 0:
            #log.debug('zone_type=%s sline=%s' %(zone_type, sline))


        if not len(sline) == len(xyz_results[inode, :]):
            msg = 'sline=%s xyz_results[inode, :]=%s' % (sline, xyz_results[inode, :])
            raise RuntimeError(msg)

        try:
            xyz_results[inode, :] = sline
        except ValueError:
            msg = 'i=%s line=%r\n' % (inode, line)
            msg += 'sline = %s' % str(sline)
            print(msg)
            raise
        iline, line, sline = get_next_sline(lines, iline)
        #log.debug(sline)
    log.debug('end of POINT')
    return iline, line, sline

def read_block(lines: list[str], iline: int,
               xyz_results: np.ndarray,
               zone_type: str, line: str,
               sline: list[str],
               nnodes: int, nvars: int,
               log: SimpleLogger) -> tuple[int, str, list[str]]: # (iline, line, sline):
    """
    BLOCK format is similar to PLOT3D in that you read all the X values before the Ys,
    Zs, and results.  The alternative format is POINT, which reads them on a per node
    basis.
    """
    log.debug('start of BLOCK')
    #print('nnodes =', nnodes)
    #print('nvars =', nvars)
    ndata = nnodes * nvars
    #print('ndata =', ndata)
    results = []

    while len(results) < ndata:
        sline = split_line(line)
        results += sline
        #print('block:', iline, sline, len(results))
        if len(sline) == 0:
            raise
        iline, line, sline = get_next_sline(lines, iline)
        #log.debug(sline)
    #print(len(results))
    assert len(results) == ndata, 'len(results)=%s expected=%s' % (len(results), ndata)
    log.debug('end of BLOCK')

    #TODO: save results
    raise RuntimeError('not done...save results')
    return iline, line, sline
    #return iline


def get_next_line(lines: list[str], iline: int) -> tuple[int, Optional[str]]:
    """Read the next line from the file.  Handles comments."""
    try:
        line = lines[iline].strip()
    except IndexError:
        line = None
        return iline, line

    line = cast(str, line)
    iline += 1
    igap = 0
    ngap_max = 10
    while len(line) == 0 or line[0] == '#':
        try:
            line = lines[iline].strip()
        except IndexError:
            line = None
            return iline, line
        iline += 1
        if igap > ngap_max:
            break
        igap += 1
    return iline, line

def get_next_sline(lines: list[str], iline: int,
                   ) -> tuple[int, Optional[str], Optional[list[str]]]:
    """Read the next split line from the file.  Handles comments."""
    iline, line = get_next_line(lines, iline)
    if line is None:
        return iline, None, None
    sline = split_line(line)
    return iline, line, sline

def get_next_nsline(lines: list[str],
                    iline: int,
                    sline: list[str],
                    nvars: int) -> tuple[int, list[str]]:
    #print(iline, sline)
    while len(sline) != nvars:  # long line was split
        #print(sline, nvars)
        iline, line, slinei = get_next_sline(lines, iline)
        #print(iline, line, slinei, nvars)
        assert len(slinei) > 0, slinei
        sline += slinei
        #print(sline, '\n')
        #iline += 1
    assert len(sline) == nvars, 'iline=%i sline=%s nvars=%s' % (iline, sline, nvars)
    return iline, sline

def split_line(line: str) -> list[str]:
    """splits a comma or space separated line"""
    if ',' in line:
        line2 = line.replace(',', ' ')
        sline = line2.split()
    else:
        sline = line.split()
    return sline

