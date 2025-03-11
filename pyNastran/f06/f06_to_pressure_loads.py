from collections import defaultdict
from typing import Optional

from cpylog import SimpleLogger
from pyNastran.bdf.bdf import read_bdf, print_card_8
from pyNastran.f06.parse_trim import read_f06_trim


def f06_to_pressure_loads(f06_filename: str,
                          subpanel_caero_filename: str,
                          loads_filename: str,
                          log: Optional[SimpleLogger]=None,
                          nlines_max: int=1_000_000,
                          debug: bool=False):
    caero_model = read_bdf(subpanel_caero_filename, xref=False, validate=False)

    nid_to_eid_map = defaultdict(list)
    for eid, elem in caero_model.elements.items():
        for nid in elem.nodes:
            nid_to_eid_map[nid].append(eid)

    trim_results = read_f06_trim(
        f06_filename, nlines_max=nlines_max,
        log=log, debug=debug)['trim_results']

    metadata = trim_results.metadata
    #print('trim_results.aero_pressure', trim_results.aero_pressure)

    with open(loads_filename, 'w') as loads_file:
        for subcase, datai in trim_results.aero_pressure.items():
            element_pressures = defaultdict(list)
            metadatai = metadata[subcase]
            subtitle = metadatai.get('subtitle', '')
            mach = metadatai['mach']
            q = metadatai['q']
            cref = metadatai['cref']
            bref = metadatai['bref']
            sref = metadatai['sref']
            nids, cp_pressure = datai
            #                    AERODYNAMIC PRES.       AERODYNAMIC
            # GRID   LABEL          COEFFICIENTS           PRESSURES             ELEMENT
            # 156     LS            6.035214E-01         9.052821E-01            900014
            for nid, (cp, p) in zip(nids, cp_pressure):
                #print(nid, cp, q)
                eids = nid_to_eid_map[nid]
                for eid in eids:
                    element_pressures[eid].append(cp)
            comment = f'$ subtitle={subtitle!r}\n'
            comment += f'$ mach={mach:g} q={q:g}\n'
            comment += f'$ bref={cref:g} bref{bref:g} sref={sref:g}\n'
            loads_file.write(comment)
            for eid, cps in element_pressures.items():
                cp = np.mean(cps)
                card = ['PLOAD2', subcase, eid, cp]
                loads_file.write(print_card_8(card))

    #print(out)
    #tables = out['tables']
