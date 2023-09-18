from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

#from .utils import get_group_name
if TYPE_CHECKING:
    from pyNastran.dev.bdf_vectorized3.bdf import BDF
    from tables import Group


def load_h5_material(model: BDF, input_group: Group):
    for h5_element in input_group._f_iter_nodes():
        name = h5_element.name
        #print(f'loading {name}')
        data = h5_element.read()
        material_id = data['MID']
        nmaterials = len(material_id)

        if name == 'MAT1':
            mat = model.mat1
            #dtype([('MID', '<i8'), ('E', '<f8'), ('G', '<f8'), ('NU', '<f8'),
                   #('RHO', '<f8'), ('A', '<f8'), ('TREF', '<f8'), ('GE', '<f8'),
                   #('ST', '<f8'), ('SC', '<f8'), ('SS', '<f8'), ('MCSID', '<i8'),
                   #('DOMAIN_ID', '<i8')])
            mat.material_id = material_id
            mat.E = data['E']
            mat.G = data['G']
            mat.nu = data['NU']
            mat.rho = data['RHO']
            mat.alpha = data['A']
            mat.tref = data['TREF']
            mat.St = data['ST']
            mat.Sc = data['SC']
            mat.Ss = data['SS']
            mat.mcsid = data['MCSID']
        elif name == 'MAT2':
            mat = model.mat2
            # dtype([('MID', '<i8'), ('G11', '<f8'), ('G12', '<f8'), ('G13', '<f8'),
            #        ('G22', '<f8'), ('G23', '<f8'), ('G33', '<f8'), ('RHO', '<f8'),
            #        ('A1', '<f8'), ('A2', '<f8'), ('A12', '<f8'), ('TREF', '<f8'),
            #        ('GE', '<f8'), ('ST', '<f8'), ('SC', '<f8'), ('SS', '<f8'),
            #        ('MCSID', '<i8'), ('DOMAIN_ID', '<i8')])
            G11 = data['G11']
            G12 = data['G12']
            G13 = data['G13']
            G22 = data['G22']
            G23 = data['G23']
            G33 = data['G33']
            rho = data['RHO']
            alpha = np.hstack([
                data['A1'].reshape(nmaterials, 1),
                data['A2'].reshape(nmaterials, 1),
                data['A12'].reshape(nmaterials, 1),
            ])
            tref = data['TREF']
            St = data['ST']
            Sc = data['SC']
            Ss = data['SS']
            ge = data['GE']
            mcsid = data['MCSID']
            mat._save(material_id, G11, G12, G13, G22, G23, G33, rho,
                      alpha, tref, ge, Ss, St, Sc, mcsid, ge_matrix=None)
        elif name == 'MAT8':
            mat = model.mat8
            # dtype([('MID', '<i8'), ('E1', '<f8'), ('E2', '<f8'), ('NU12', '<f8'),
            #        ('G12', '<f8'), ('G1Z', '<f8'), ('G2Z', '<f8'), ('RHO', '<f8'),
            #        ('A1', '<f8'), ('A2', '<f8'), ('TREF', '<f8'), ('XT', '<f8'),
            #        ('XC', '<f8'), ('YT', '<f8'), ('YC', '<f8'), ('S', '<f8'),
            #        ('GE', '<f8'), ('F12', '<f8'), ('STRN', '<f8'), ('DOMAIN_ID', '<i8')])

            E1 = data['E1']
            E2 = data['E2']
            nu12 = data['NU12']
            G12 = data['G12']
            G1z = data['G1Z']
            G2z = data['G2Z']
            rho = data['RHO']
            alpha = np.hstack([
                data['A1'].reshape(nmaterials, 1),
                data['A2'].reshape(nmaterials, 1),
            ])
            tref = data['TREF']
            Xt = data['XT']
            Xc = data['XC']
            Yt = data['YT']
            Yc = data['YC']
            S = data['S']
            ge = data['GE']
            f12 = data['F12']
            strn = data['STRN']

            mat._save(material_id, E1, E2, G12, G1z, G2z, nu12,
                      rho, alpha, tref, ge, Xt, Xc, Yt, Yc, S, f12, strn)
        else:
            raise NotImplementedError(name)
        material_id = material_id
        mat.n = nmaterials
        mat.write()
        x = 1
    x = 2
