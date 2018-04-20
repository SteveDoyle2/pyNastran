from __future__ import print_function, absolute_import

import numpy as np
from numpy import cos, sin
from pyNastran.op2.data_in_material_coord import transf_Mohr


from ..result_table import ResultTable


def rotate_shell_loading_force_stress(shell_loading):
    nx, ny, nxy, mx, my, mxy, vx, vy, theta_rad = shell_loading
    nx[:], ny[:], nxy[:] = transf_Mohr(nx, ny, nxy, theta_rad)
    mx[:], my[:], mxy[:] = transf_Mohr(mx, my, mxy, theta_rad)

    c = cos(theta_rad)
    s = sin(theta_rad)

    vx[:], vy[:] = c * vx + s * vy, -s * vx + c * vy
    
    
def rotate_shell_loading_force_stress_complex(shell_loading):
    shell_loading_real = shell_loading[:8] + shell_loading[-1]
    rotate_shell_loading_force_stress(shell_loading_real)
    
    shell_loading_complex = shell_loading[8:]
    rotate_shell_loading_force_stress(shell_loading_complex)
    
    
def rotate_shell_loading_strain(strains):
    ex, ey, exy, theta_rad = strains
    exy /= 2.
    ex[:], ey[:], exy[:] = transf_Mohr(ex, ey, exy, theta_rad)
    exy *= 2.


def rotate_shell_loading_strain_complex(strains):
    shell_loading_real = strains[:3] + strains[-1]
    rotate_shell_loading_strain(shell_loading_real)

    shell_loading_complex = strains[3:]
    rotate_shell_loading_strain(shell_loading_complex)
    

class ShellElementForceStressResultTable(object):
    def apply_options(self, data):
        if 'MATERIAL' in data.options:
            _data = data.data
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in _data['EID']])
            shell_loading = [_data[col] for col in self.result_data_cols] + [-theta_rad]
            rotate_shell_loading_force_stress(shell_loading)

    def search(self, data_ids, domains=(), material_sys=False):
        result = ResultTable.search(self, data_ids, domains)

        if material_sys:
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in result['EID'].values])
            shell_loading = [result[col].values for col in self.result_data_cols] + [theta_rad]
            rotate_shell_loading_force_stress(shell_loading)

        return result
    
    
class ShellElementForceStressResultTableComplex(object):
    def apply_options(self, data):
        if 'MATERIAL' in data.options:
            _data = data.data
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in _data['EID']])
            shell_loading = [_data[col] for col in self.result_data_cols] + [-theta_rad]
            rotate_shell_loading_force_stress_complex(shell_loading)

    def search(self, data_ids, domains=(), material_sys=False):
        result = ResultTable.search(self, data_ids, domains)

        if material_sys:
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in result['EID'].values])
            shell_loading = [result[col].values for col in self.result_data_cols] + [theta_rad]
            rotate_shell_loading_force_stress_complex(shell_loading)

        return result


class ShellElementStrainResultTable(object):
    def apply_options(self, data):
        if 'MATERIAL' in data.options:
            _data = data.data
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in _data['EID']])
            strains = [_data[col] for col in self.result_data_cols] + [-theta_rad]
            rotate_shell_loading_strain(strains)

    def search(self, data_ids, domains=(), material_sys=False):
        result = ResultTable.search(self, data_ids, domains)

        if material_sys:
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in result['EID'].values])
            strains = [result[col].values for col in self.result_data_cols] + [theta_rad]
            rotate_shell_loading_strain(strains)

        return result
    
    
class ShellElementStrainResultTableComplex(object):
    def apply_options(self, data):
        if 'MATERIAL' in data.options:
            _data = data.data
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in _data['EID']])
            strains = [_data[col] for col in self.result_data_cols] + [-theta_rad]
            rotate_shell_loading_strain_complex(strains)

    def search(self, data_ids, domains=(), material_sys=False):
        result = ResultTable.search(self, data_ids, domains)

        if material_sys:
            theta_rad = self._h5n.nastran.input.element.get_shell_element_info_dict()['THETA_RAD']
            theta_rad = np.array([theta_rad[eid] for eid in result['EID'].values])
            strains = [result[col].values for col in self.result_data_cols] + [theta_rad]
            rotate_shell_loading_strain_complex(strains)

        return result
