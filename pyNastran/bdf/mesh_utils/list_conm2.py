from __future__ import annotations
import warnings
from typing import Optional
from cpylog import SimpleLogger
import numpy as np

from pyNastran.utils import PathLike
from .cmd_line.utils_bdf import read_lax_bdf


def list_conm2(bdf_filename: PathLike,
               mass_scale: float=1.0,
               log: Optional[SimpleLogger]=None):
    cards_to_cards = [
        'GRID', 'CONM2',
        'CORD1R', 'CORD1S', 'CORD1C',
        'CORD2R', 'CORD2S', 'CORD2C',
    ]
    model = read_lax_bdf(
        bdf_filename, punch=False, validate=True,
        xref=False, is_strict_card_parser=False,
        cards_to_read=cards_to_cards, log=log)

    nids_list = []
    eids_list = []
    mass_list = []
    for eid, elem in model.masses.items():
        if elem.type == 'CONM2':
            eids_list.append(eid)
            nids_list.append(elem.nid)
            mass_list.append(elem.mass)
        else:
            warnings.warn(f'skipping {elem.type}\n{str(elem)}')

    nids = np.array(nids_list)
    eids = np.array(eids_list)
    mass = np.array(mass_list)
    imass = np.argsort(mass)
    for eid, nid in zip(eids[imass], nids[imass]):
        node = model.nodes[nid]
        elem = model.masses[eid]
        elem.mass *= mass_scale
        elem.I *= mass_scale
        print(node)
        print(elem)
    #level = 'debug' if not quiet else 'warning'
    #log = SimpleLogger(level=level, encoding='utf-8')
