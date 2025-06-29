from collections import defaultdict
from typing import Optional
import numpy as np

from cpylog import SimpleLogger
from pyNastran.bdf.bdf import read_bdf, print_card_8
#from pyNastran.f06.f06_tables.trim import AeroPressure
from pyNastran.f06.parse_trim import read_f06_trim


def f06_to_pressure_loads(f06_filename: str,
                          subpanel_caero_filename: str,
                          loads_filename: str,
                          nid_csv_filename: str='',
                          eid_csv_filename: str='',
                          log: Optional[SimpleLogger]=None,
                          nlines_max: int=1_000_000,
                          debug: bool=False) -> None:
    caero_model = read_bdf(subpanel_caero_filename, log=log,
                           xref=False, validate=False, debug=debug)
    log = caero_model.log

    nid_to_eid_map = defaultdict(list)
    for eid, elem in caero_model.elements.items():
        for nid in elem.nodes:
            nid_to_eid_map[nid].append(eid)

    trim_results = read_f06_trim(
        f06_filename, nlines_max=nlines_max,
        log=log, debug=debug)['trim_results']

    metadata = trim_results.metadata
    #print('trim_results.aero_pressure', trim_results.aero_pressure)

    element_pressures = {}
    for subcase, apress in trim_results.aero_pressure.items():
        element_pressure = apress.get_element_pressure(nid_to_eid_map)
        element_pressures[subcase] = element_pressure

    import sys
    if loads_filename is not None:
        with open(loads_filename, 'w') as loads_file:
            for subcase, element_pressure in element_pressures.items():
                if 0:  # pragma: no cover
                    metadatai = metadata[subcase]
                    subtitle = metadatai.get('subtitle', '')
                    mach = metadatai['mach']
                    q = metadatai['q']
                    cref = metadatai['cref']
                    bref = metadatai['bref']
                    sref = metadatai['sref']
                else:
                    subtitle = element_pressure.subtitle
                    mach = element_pressure.mach
                    q = element_pressure.q
                    cref = element_pressure.chord
                    bref = element_pressure.span
                    sref = element_pressure.sref

                comment = f'$ subtitle={subtitle!r}\n'
                comment += f'$ mach={mach:g} q={q:g}\n'
                comment += f'$ bref={cref:g} bref{bref:g} sref={sref:g}\n'
                loads_file.write(comment)
                for eid, cps in element_pressure.items():
                    cp = np.mean(cps)
                    card = ['PLOAD2', subcase, eid, cp]
                    loads_file.write(print_card_8(card))
        log.info(f'finished writing {loads_filename}')

    if nid_csv_filename:
        node_line0 = '# Nid,'
        for subcase, element_pressure in element_pressures.items():
            node_line0 += f'CpSubcase{subcase:d}(f),'
            eids = list(element_pressure)
        node_line0 += '\n'

        nodes = apress.nodes
        nnodes = len(nodes)
        nsubcases = len(element_pressures)
        isubcase = 0
        node_cp_array = np.zeros((nnodes, nsubcases))
        for subcase, apress in trim_results.aero_pressure.items():
            node_cp_array[:, isubcase] = apress.cp
            isubcase += 1

        with open(nid_csv_filename, 'w') as csv_file:
            csv_file.write(node_line0)
            for nid, cp_arrayi in zip(nodes, node_cp_array):
                data = [nid] + list(cp_arrayi)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {nid_csv_filename}')

    if eid_csv_filename:
        line0 = '# Eid,'
        cps_list = []
        for subcase, element_pressure in element_pressures.items():
            line0 += f'CpSubcase{subcase:d}(f),'
            eids = list(element_pressure)
            neids = len(eids)
            cp_list = []
            for eid, cps in element_pressure.items():
                cpi = np.mean(cps)
                cp_list.append(cpi)
            cps_list.append(cp_list)
            line0 = line0.rstrip(',') + '\n'
        cp_array = np.column_stack(cps_list)
        assert cp_array.shape == (neids, nsubcases), f'actual_shape={cp_array.shape} expected=({neids},{nsubcases})'

        with open(eid_csv_filename, 'w') as csv_file:
            csv_file.write(line0)
            for eid, cp_arrayi in zip(eids, cp_array):
                data = [eid] + list(cp_arrayi)
                strs = ','.join(map(str, data))
                csv_file.write(strs + '\n')
        log.info(f'finished writing {eid_csv_filename}')
    #print(out)
    #tables = out['tables']
