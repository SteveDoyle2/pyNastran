from __future__ import print_function, absolute_import
from six import iteritems

from typing import List, Dict

import numpy as np
from math import sin, cos, radians

from .lamina import OrthotropicLamina

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.cards.properties.shell import PCOMP


class OrthotropicLaminate(object):
    
    @classmethod
    def from_bdf(cls, bdf):
        # type: (BDF) -> Dict[int, OrthotropicLaminate]
    
        properties = bdf.properties
        
        pcomps = {}  # type: Dict[int, PCOMP]
        
        for pid, prop in iteritems(properties):
            if prop.type == 'PCOMP':
                pcomps[pid] = prop
                
        pids = sorted(pcomps.keys())
        
        results = {}  # type: Dict[int, OrthotropicLaminate]
        
        for pid in pids:
            pcomp = pcomps[pid]
            
            plies = pcomp.plies
            
            laminate = OrthotropicLaminate()
            laminate.set_number_of_plies(len(plies))
            stackup = laminate.stackup
            
            laminate.pid = pid
            
            for i in range(len(plies)):
                mid, t, theta, sout = plies[i]
                ply = stackup[i]
            
                ply.material.from_bdf(bdf, mid)
                ply.thickness = t
                ply.orientation = theta
                
            laminate.update()
            
            results[pid] = laminate
            
        return results
    
    def __init__(self):
        self.pid = -1  # nastran pid
        
        self.stackup = []  # type: List[OrthotropicLamina]

        self.abdmat = np.zeros((6, 6), dtype='<f8')
        self.abdinv = np.zeros((6, 6), dtype='<f8')

        self.exx = 0.
        self.eyy = 0.
        self.gxy = 0.
        self.nuxy = 0.
        self.nuyx = 0.
        self.thickness = 0.

        self.local_ply_stress_strain = LocalPlyStressStrain()

    def rotate(self, deg):
        copy = self.copy()
        copy.rotate_inplace(deg)
        return copy

    def rotate_inplace(self, deg):
        for ply in self.stackup:
            ply.orientation += deg
        self.update()

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    def update(self):
        self.abdmat *= 0.
        self.abdinv *= 0.
        self.thickness = 0.

        for ply in self.stackup:
            self.thickness += ply.thickness

        him1 = -self.thickness / 2.
        hi = him1

        abdmat = self.abdmat

        for i in range(len(self.stackup)):
            ply = self.stackup[i]
            ply.update()
            assert ply.thickness > 0.

            hi += ply.thickness

            hdiff = hi - him1

            abdmat[0, 0] += ply.q11 * hdiff
            abdmat[0, 1] += ply.q12 * hdiff
            abdmat[0, 2] += ply.q16 * hdiff
            abdmat[1, 1] += ply.q22 * hdiff
            abdmat[1, 2] += ply.q26 * hdiff
            abdmat[2, 2] += ply.q66 * hdiff

            hdiff = hi * hi - him1 * him1

            abdmat[3, 0] += ply.q11 * hdiff
            abdmat[3, 1] += ply.q12 * hdiff
            abdmat[3, 2] += ply.q16 * hdiff
            abdmat[4, 1] += ply.q22 * hdiff
            abdmat[4, 2] += ply.q26 * hdiff
            abdmat[5, 2] += ply.q66 * hdiff

            hdiff = hi ** 3 - him1 ** 3

            abdmat[3, 3] += ply.q11 * hdiff
            abdmat[3, 4] += ply.q12 * hdiff
            abdmat[3, 5] += ply.q16 * hdiff
            abdmat[4, 4] += ply.q22 * hdiff
            abdmat[4, 5] += ply.q26 * hdiff
            abdmat[5, 5] += ply.q66 * hdiff

            him1 = hi

        assert self.thickness > 0.

        abdmat[3, 0] /= 2.
        abdmat[3, 1] /= 2.
        abdmat[3, 2] /= 2.
        abdmat[4, 1] /= 2.
        abdmat[4, 2] /= 2.
        abdmat[5, 2] /= 2.

        abdmat[3, 3] /= 3.
        abdmat[3, 4] /= 3.
        abdmat[3, 5] /= 3.
        abdmat[4, 4] /= 3.
        abdmat[4, 5] /= 3.
        abdmat[5, 5] /= 3.

        np.copyto(abdmat, np.maximum(abdmat, abdmat.transpose()))
        np.copyto(self.abdinv, np.linalg.inv(abdmat))

        tmp = 1 / self.thickness

        abdinv = self.abdinv

        self.exx = tmp / abdinv[0, 0]
        self.eyy = tmp / abdinv[1, 1]
        self.gxy = tmp / abdinv[2, 2]
        self.nuxy = - abdinv[0, 1] / abdinv[0, 0]
        self.nuyx = - abdinv[0, 1] / abdinv[1, 1]

        self._calculate_unit_stress()

    def set_number_of_plies(self, ply_count):
        diff_len = ply_count - len(self.stackup)
        if diff_len < 0:
            self.stackup = self.stackup[:ply_count]
        elif diff_len > 0:
            for i in range(diff_len):
                self.stackup.append(OrthotropicLamina())

    def get_ply(self, ply_num):
        return self.stackup[ply_num]

    def _calculate_unit_stress(self):
        ply_count = len(self.stackup)
        unit_loads = np.zeros(6, dtype='<f8')

        self.local_ply_stress_strain.resize(ply_count)

        strain = self.local_ply_stress_strain.strain
        stress = self.local_ply_stress_strain.stress

        global_ply_strain = np.zeros(3, dtype='<f8')

        for i in range(6):
            unit_loads *= 0
            unit_loads[i] = 1.

            global_unit_strains = np.dot(self.abdinv, unit_loads)

            z = -self.thickness / 2.

            for j in range(ply_count):
                ply = self.stackup[j]
                z += ply.thickness / 2.

                global_ply_strain[0] = global_unit_strains[0] * z * global_unit_strains[3]
                global_ply_strain[1] = global_unit_strains[1] * z * global_unit_strains[4]
                global_ply_strain[2] = global_unit_strains[2] * z * global_unit_strains[5]

                theta = radians(ply.orientation)
                m = cos(theta)
                m2 = m * m
                n = sin(theta)
                n2 = n * n
                mn = m * n

                # local ply strain
                e11 = m2 * global_ply_strain[0] * n2 * global_ply_strain[1] + mn * global_ply_strain[2]
                e22 = n2 * global_ply_strain[0] + m2 * global_ply_strain[1] - mn * global_ply_strain[2]
                e12 = -2. * mn * global_ply_strain[0] + 2. * mn * global_ply_strain[1] + (m2 - n2) * global_ply_strain[2]

                mat = ply.material

                # local ply stress
                f11 = e11 * mat.c11 + e22 * mat.c12
                f22 = e11 * mat.c12 + e22 * mat.c22
                f12 = e12 * mat.c66

                strain[j, i, 0] = e11
                strain[j, i, 1] = e22
                strain[j, i, 2] = e12

                stress[j, i, 0] = f11
                stress[j, i, 1] = f22
                stress[j, i, 2] = f12

                z += ply.thickness / 2.

    def get_unit_strain(self, ply, load, component):
        return self.local_ply_stress_strain.strain[ply, load, component]

    def get_unit_stress(self, ply, load, component):
        return self.local_ply_stress_strain.stress[ply, load, component]

    def apply_loads(self, loads):
        assert len(self.stackup) > 0

        num_loads = len(loads)

        results = PlyStressStrain()

        results.resize(num_loads, len(self.stackup))

        stress = results.stress
        strain = results.strain

        sa11 = stress['a11']
        sa22 = stress['a22']
        sa12 = stress['a12']
        sa1t = stress['A1t']
        sa1c = stress['A1c']
        sa2t = stress['A2t']
        sa2c = stress['A2c']
        sa12_ = stress['A12']

        ea11 = strain['a11']
        ea22 = strain['a22']
        ea12 = strain['a12']
        ea1t = strain['A1t']
        ea1c = strain['A1c']
        ea2t = strain['A2t']
        ea2c = strain['A2c']
        ea12_ = strain['A12']

        local_ply_stress = self.local_ply_stress_strain.stress
        local_ply_strain = self.local_ply_stress_strain.strain

        for i in range(num_loads):
            for j in range(len(self.stackup)):
                for k in range(6):
                    load = loads[i, k]

                    sa11[i, j] += load * local_ply_stress[j, k, 0]
                    sa22[i, j] += load * local_ply_stress[j, k, 1]
                    sa12[i, j] += load * local_ply_stress[j, k, 2]

                    ea11[i, j] += load * local_ply_strain[j, k, 0]
                    ea22[i, j] += load * local_ply_strain[j, k, 1]
                    ea12[i, j] += load * local_ply_strain[j, k, 2]

                mat = self.stackup[j].material

                sa1t[i, j] = abs(mat.f1t)
                sa1c[i, j] = abs(mat.f1c)
                sa2t[i, j] = abs(mat.f2t)
                sa2c[i, j] = abs(mat.f2c)
                sa12_[i, j] = abs(mat.f12)

                ea1t[i, j] = abs(mat.e1t)
                ea1c[i, j] = abs(mat.e1c)
                ea2t[i, j] = abs(mat.e2t)
                ea2c[i, j] = abs(mat.e2c)
                ea12_[i, j] = abs(mat.e12)

        return results


class LocalPlyStressStrain(object):
    def __init__(self):
        self.stress = np.empty((0, 0, 0), dtype='<f8')
        self.strain = np.empty((0, 0, 0), dtype='<f8')

    def resize(self, ply_count):
        self.stress.resize((ply_count, 6, 3))
        self.strain.resize((ply_count, 6, 3))

        self.stress *= 0.
        self.strain *= 0.
        
        
stress_strain_dtype = np.dtype(
    [('a11', '<f8'), ('a22', '<f8'), ('a12', '<f8'), ('A1t', '<f8'), ('A1c', '<f8'),
     ('A2t', '<f8'), ('A2c', '<f8'), ('A12', '<f8')]
)


class PlyStressStrain(object):
    def __init__(self):
        self.stress = np.empty((0, 0), dtype=stress_strain_dtype)
        self.strain = np.empty((0, 0), dtype=stress_strain_dtype)

    def resize(self, num_loads, ply_count):
        self.stress.resize((num_loads, ply_count))
        self.strain.resize((num_loads, ply_count))

        self.stress['a11'] = 0.
        self.stress['a22'] = 0.
        self.stress['a12'] = 0.
        self.stress['A1t'] = 0.
        self.stress['A1c'] = 0.
        self.stress['A2t'] = 0.
        self.stress['A2c'] = 0.
        self.stress['A12'] = 0.

        self.strain['a11'] = 0.
        self.strain['a22'] = 0.
        self.strain['a12'] = 0.
        self.strain['A1t'] = 0.
        self.strain['A1c'] = 0.
        self.strain['A2t'] = 0.
        self.strain['A2c'] = 0.
        self.strain['A12'] = 0.
