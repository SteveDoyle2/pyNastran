from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from .utils import get_group_name
from tables import Group, Node, File, Table
from pyNastran.dev.bdf_vectorized3.bdf_interface.h5_pytables.utils import get_group_name
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.dev.bdf_vectorized3.bdf import BDF


def load_h5_dynamic(model: BDF, input_group: Group):
    log = model.log
    for h5_element in input_group._f_iter_nodes():
        if isinstance(h5_element, Group):
            #name = h5_node.name
            ## TODO: janky way to get the name
            name = get_group_name(h5_element)
            log.info(f'loading {name}')
            
            # groups
            if name == 'SPC1':
                assert len(h5_element.dtype) == 8, h5_element.dtype
                grids = h5_element['SPC1_G']
                thrus = h5_element['SPC1_THRU'].read()
                log.info(f'SPC1 grids={grids} thrus={thrus}')
                #for h5_elementi in h5_element._f_iter_nodes():
                #    print(f'loading {h5_elementi}')

            # tables
            elif name == 'TEMP':
                # [('SID', '<i8'), ('G', '<i8'), ('T', '<f8'), ('DOMAIN_ID', '<i8')]
                assert len(h5_element.dtype) == 4, h5_element.dtype
                table = h5_element.read()
                #log.info(f'{name} table={table}')
                for sid, nid, temp, domain_id in table:
                    model.add_temp(sid, {nid: temp})
            elif name == 'TEMPD':
                # [('SID', '<i8'), ('T', '<f8'), ('DOMAIN_ID', '<i8')]
                assert len(h5_element.dtype) == 3, h5_element.dtype
                table = h5_element.read()
                #log.info(f'{name} table={table}')
                for sid, temp, domain_id in table:
                    model.add_tempd(sid, temp)
            else:
                # group
                log.info(f'loading {name}')
                for h5_elementi in h5_element._f_iter_nodes():
                    log.debug(f'  {h5_elementi}')
        elif isinstance(h5_element, Table):
            name = get_group_name(h5_element)
            name = get_group_name(h5_element)
            log.info(f'loading {name}')
            out = h5_element.read()
            log.info(out)
            log.info(h5_element)
        else:
            raise RuntimeError(h5_element)
        continue
        data = h5_element.read()
        print(data)
        name = h5_element.name
        print(f'loading {name}')
        #material_id = data['MID']
        #nmaterials = len(material_id)

        if name == 'SPC1':
            mat = model.mat1
            _load_mat1(mat, data, material_id, nmaterials)
        else:  # pragma: no cover
            raise NotImplementedError(name)
        material_id = material_id
        mat.n = nmaterials
        try:
            mat.write()
        except:  # pragma: no cover
            print(mat.get_stats())
            raise
        x = 1
    x = 2


def _load_mat1(mat: MAT1, data, material_id, nmaterials: int):
    # dtype([('MID', '<i8'), ('E', '<f8'), ('G', '<f8'), ('NU', '<f8'),
    # ('RHO', '<f8'), ('A', '<f8'), ('TREF', '<f8'), ('GE', '<f8'),
    # ('ST', '<f8'), ('SC', '<f8'), ('SS', '<f8'), ('MCSID', '<i8'),
    # ('DOMAIN_ID', '<i8')])
    E = data['E']
    G = data['G']
    nu = data['NU']
    rho = data['RHO']
    alpha = data['A']
    tref = data['TREF']
    ge = data['GE']
    St = data['ST']
    Sc = data['SC']
    Ss = data['SS']
    mcsid = data['MCSID']
    mat._save(material_id, E, G, nu, rho, alpha, tref, ge,
              Ss, St, Sc, mcsid)
