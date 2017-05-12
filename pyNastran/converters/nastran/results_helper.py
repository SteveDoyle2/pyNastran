from __future__ import print_function
from copy import deepcopy
from collections import defaultdict
import traceback

from six import iteritems
import numpy as np
from numpy.linalg import norm

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.nastran.geometry_helper import NastranGuiAttributes
from pyNastran.converters.nastran.displacements import (
    DisplacementResults, ForceTableResults) #, TransientElementResults

class NastranGuiResults(NastranGuiAttributes):
    """
    Defines OP2 specific methods NastranIO
    """
    def __init__(self):
        super(NastranGuiResults, self).__init__()

    def _fill_gpforces(self, model):
        pass
        #[model.grid_point_forces, 'GridPointForces'],  # TODO: this is not really an OUG table

    def _fill_op2_oug_oqg(self, cases, model, key, icase,
                          form_dict, header_dict):
        """
        loads the nodal dispalcements/velocity/acceleration/eigenvector/spc/mpc forces
        """
        new_cases = True
        nnodes = self.nNodes
        displacement_like = [
            # slot, name, deflects

            # TODO: what is a velocity/acceleration?
            #       is it a fringe, displacement, force?
            (model.displacements, 'Displacement', True),
            (model.velocities, 'Velocity', False),
            (model.accelerations, 'Acceleration', False),
            (model.eigenvectors, 'Eigenvectors', True),
            (model.spc_forces, 'SPC Forces', False),
            (model.mpc_forces, 'MPC Forces', False),

            (model.load_vectors, 'LoadVectors', False),
            (model.applied_loads, 'AppliedLoads', False),
            (model.force_vectors, 'ForceVectors', False),
        ]
        temperature_like = [
            (model.temperatures, 'Temperature'),
        ]
        nids = self.node_ids
        for (result, name, deflects) in displacement_like:
            if key not in result:
                continue

            title1 = name + ' T_XYZ'
            #title2 = name + ' R_XYZ'

            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                print('str(%s) has no data...' % case.__class.__name__)
                continue
            #if not case.is_real():
                #print('complex results not supported...')
                #continue
            # transient
            if case.nonlinear_factor is not None:
                #code_name = case.data_code['name']
                has_cycle = hasattr(case, 'mode_cycle')
            else:
                has_cycle = False
                code_name = None
            assert case.is_sort1(), case.is_sort1()

            itime0 = 0
            t1 = case.data[itime0, :, 0]
            ndata = t1.shape[0]
            if nnodes != ndata:
                #print('nnodes=%s ndata=%s' % (nnodes, ndata))
                nidsi = case.node_gridtype[:, 0]
                assert len(nidsi) == nnodes
                j = np.searchsorted(nids, nidsi)  # searching for nidsi

                try:
                    if not np.allclose(nids[j], nidsi):
                        msg = 'nids[j]=%s nidsi=%s' % (nids[j], nidsi)
                        raise RuntimeError(msg)
                except IndexError:
                    msg = 'node_ids = %s\n' % list(nids)
                    msg += 'nidsi in disp = %s\n' % list(nidsi)
                    raise IndexError(msg)

            # (itime, nnodes, xyz)
            # (901, 6673, 3)
            t123 = case.data[:, :, :3]
            if nnodes != ndata:
                t123i = np.zeros((nnodes, 3), dtype='float32')
                t123i[j, :] = t123
                t123 = t123i

            # (itime, nnodes, xyz)
            # tnorm (901, 3)
            tnorm = norm(t123, axis=1)   # I think this is wrong...
            assert len(tnorm) == t123.shape[0]
            ntimes = case.ntimes
            titles = []
            scales = []
            headers = []
            #if deflects:
            if deflects:
                nastran_res = DisplacementResults(subcase_idi, titles, headers,
                                                  self.xyz_cid0, t123, tnorm,
                                                  scales, #deflects=deflects,
                                                  uname='NastranResult')

                dmax = []
                for itime in range(ntimes):
                    dt = case._times[itime]

                    if name == 'Displacement':
                        # (6673, )
                        normiii = np.linalg.norm(t123[itime, :, :], axis=1)
                        #print(normiii.shape)
                        #print('Displacement; itime=%s time=%s tnorm=%s' % (itime, dt, normiii.max()))
                        dmax.append(normiii.max())
                    # mode = 2; freq = 75.9575 Hz
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header

                    tnorm_abs_max = tnorm.max()
                    #if tnorm_abs_max == 0.0:
                        #scale = self.displacement_scale_factor
                    #else:
                        #scale = self.displacement_scale_factor / tnorm_abs_max

                    scale = self.dim_max
                    if tnorm_abs_max > 0.0:
                        scale = self.dim_max / tnorm_abs_max * 0.25
                    scales.append(scale)
                    titles.append(title1)
                    headers.append(header)
                    cases[icase] = (nastran_res, (itime, title1))  # do I keep this???
                    formii = (title1, icase, [])
                    form_dict[(key, itime)].append(formii)
                    icase += 1

                if name == 'Displacement':
                    # Displacement; itime=361 time=3.61 tnorm=1.46723
                    print('dmax = ', max(dmax))
                nastran_res.save_defaults()
            else:
                nastran_res = ForceTableResults(subcase_idi, titles, headers,
                                                t123, tnorm,
                                                scales, #deflects=deflects,
                                                uname='NastranResult')
                for itime in range(ntimes):
                    dt = case._times[itime]
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header

                    tnorm_abs_max = tnorm.max()
                    #if tnorm_abs_max == 0.0:
                        #scale = self.displacement_scale_factor
                    #else:
                        #scale = self.displacement_scale_factor / tnorm_abs_max

                    # TODO: what to do with the scale factor?
                    #scale = self.dim_max
                    #if tnorm_abs_max > 0.0:
                        #scale = self.dim_max / tnorm_abs_max * 0.25
                    scale = 1.
                    scales.append(scale)
                    titles.append(title1)
                    headers.append(header)
                    cases[icase] = (nastran_res, (itime, title1))  # do I keep this???
                    formii = (title1, icase, [])
                    form_dict[(key, itime)].append(formii)
                    icase += 1
                nastran_res.save_defaults()

        for (result, name) in temperature_like:
            if key not in result:
                continue
            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                continue

            ntimes = case.ntimes
            for itime in range(ntimes):
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header

                loads = case.data[itime, :, :]
                nxyz = norm(loads[:, :3], axis=1)
                assert len(nxyz) == nnodes, 'len(nxyz)=%s nnodes=%s' % (
                    len(nxyz), nnodes)

                temp_res = GuiResult(subcase_idi, header=name, title=name,
                                     location='node', scalar=loads[:, 0])
                cases[icase] = (temp_res, (0, name))
                form_dict[(key, itime)].append((name, icase, []))
                icase += 1
        return icase

    def _fill_op2_force2(self, cases, model, key, icase, itime,
                         form_dict, header_dict, is_static):
        """creates the thermal loads"""
        thermal_loads = [
            # 3D
            model.chexa_thermal_load, model.ctetra_thermal_load,
            model.cpenta_thermal_load,
        ]
        eids = self.element_ids
        name = 'thermal_load'
        for result in thermal_loads:
            if key not in result:
                continue

            title = name + 'XYZ'
            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                print('str(%s) has no data...' % case.__class.__name__)
                continue
            #if not case.is_real():
                #print('complex results not supported...')
                #continue
            # transient
            if case.nonlinear_factor is not None:
                #code_name = case.data_code['name']
                has_cycle = hasattr(case, 'mode_cycle')
            else:
                has_cycle = False
                code_name = None
            assert case.is_sort1(), case.is_sort1()

            itime0 = 0
            t1 = case.data[itime0, :, 0]
            ndata = t1.shape[0]
            if nelements != ndata:
                eidsi = case.elements
                assert len(eidsi) == nelements
                j = np.searchsorted(eids, eidsi)  # searching for eidsi

                try:
                    if not np.allclose(eids[j], eidsi):
                        msg = 'nids[j]=%s eidsi=%s' % (eids[j], eidsi)
                        raise RuntimeError(msg)
                except IndexError:
                    msg = 'element_ids = %s\n' % list(eids)
                    msg += 'eidsi in force = %s\n' % list(eidsi)
                    raise IndexError(msg)

            # (itime, nelements, xyz)
            t123 = case.data[:, :, :3]
            if nelements != ndata:
                t123i = np.zeros((nelements, 3), dtype='float32')
                t123i[j, :] = t123
                t123 = t123i
            tnorm = norm(t123, axis=1)
            assert len(tnorm) == t123.shape[0]
            ntimes = case.ntimes
            titles = []
            scales = []
            headers = []
            for itime in range(ntimes):
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header

                loads = case.data[itime, :, :]
                txyz = norm(loads[:, :3], axis=1)
                rxyz = norm(loads[:, 3:6], axis=1)
                assert loads[:, :3].shape[1] == 3, loads.shape
                assert loads[:, 3:6].shape[1] == 3, loads.shape
                assert len(txyz) == nelements, 'len(txyz)=%s nnodes=%s' % (
                    len(txyz), nnodes)

                tx_res = GuiResult(subcase_idi, header=name + 'Tx', title=name + 'Tx',
                                   location='node', scalar=loads[:, 0])
                ty_res = GuiResult(subcase_idi, header=name + 'Ty', title=name + 'Ty',
                                   location='node', scalar=loads[:, 1])
                tz_res = GuiResult(subcase_idi, header=name + 'Tz', title=name + 'Tz',
                                   location='node', scalar=loads[:, 2])
                rx_res = GuiResult(subcase_idi, header=name + 'Rx', title=name + 'Rx',
                                   location='node', scalar=loads[:, 3])
                ry_res = GuiResult(subcase_idi, header=name + 'Ry', title=name + 'Ry',
                                   location='node', scalar=loads[:, 4])
                rz_res = GuiResult(subcase_idi, header=name + 'Rz', title=name + 'Rz',
                                   location='node', scalar=loads[:, 5])
                txyz_res = GuiResult(subcase_idi, header=name + 'Txyz',
                                     title=name + 'Txyz', location='node', scalar=txyz)
                rxyz_res = GuiResult(subcase_idi, header=name + 'Rxyz',
                                     title=name + 'Rxyz', location='node', scalar=rxyz)

                cases[icase] = (tx_res, (0, name + 'Tx'))
                cases[icase + 1] = (ty_res, (0, name + 'Ty'))
                cases[icase + 2] = (tz_res, (0, name + 'Tz'))
                cases[icase + 3] = (txyz_res, (0, name  + 'Txyz'))
                cases[icase + 4] = (rx_res, (0, name + 'Rx'))
                cases[icase + 5] = (ry_res, (0, name + 'Ry'))
                cases[icase + 6] = (rz_res, (0, name + 'Rz'))
                cases[icase + 7] = (rxyz_res, (0, name  + 'Rxyz'))

                form_dict[(key, itime)].append((name + 'Tx', icase, []))
                form_dict[(key, itime)].append((name + 'Ty', icase + 1, []))
                form_dict[(key, itime)].append((name + 'Tz', icase + 2, []))
                form_dict[(key, itime)].append((name + 'Txyz', icase + 3, []))
                form_dict[(key, itime)].append((name + 'Rx', icase + 4, []))
                form_dict[(key, itime)].append((name + 'Ry', icase + 5, []))
                form_dict[(key, itime)].append((name + 'Rz', icase + 6, []))
                form_dict[(key, itime)].append((name + 'Rxyz', icase + 7, []))
                icase += 7

    def _fill_op2_force(self, cases, model, key, icase, itime,
                        form_dict, header_dict, is_static):
        """creates the force plots"""
        #assert isinstance(key, int), key
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        try:
            icase = self._fill_op2_time_centroidal_force(
                cases, model, key, icase, itime,
                form_dict, header_dict, is_static)
        except IndexError as e:
            self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(e))
        return icase

    def _fill_op2_stress(self, cases, model, key, icase, itime,
                         form_dict, header_dict, is_static, is_stress=True):
        """creates the stress plots"""
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        try:
            icase = self._fill_op2_time_centroidal_stress(
                cases, model, key, icase, itime, form_dict, header_dict,
                is_static, is_stress=is_stress)
        except TypeError as e:
            self.log_error('\n' + ''.join(traceback.format_stack()))
            #traceback.print_exc(file=self.log_error)
            self.log_error(str(e))
        return icase

    def _fill_op2_strain(self, cases, model, key, icase, itime,
                         form_dict, header_dict, is_static):
        """creates the strain plots"""
        return self._fill_op2_stress(cases, model, key, icase, itime,
                                     form_dict, header_dict,
                                     is_static, is_stress=False)

    def _get_times(self, model, isubcase):
        """
        Get the times/frequencies/eigenvalues/loadsteps used on a given
        subcase
        """
        table_types = model.get_table_types()
        is_real = True
        is_data = False
        is_static = False
        times = None
        #print('isubcase =', isubcase)
        #print('table_types =', table_types)
        #print('model.eigenvectors.keys() =', model.eigenvectors.keys())
        for table_type in table_types:
            if not hasattr(model, table_type):
                print('no table_type=%s' % table_type)
                continue
            table = getattr(model, table_type)
            if len(table) == 0:
                continue
            #print(table)
            if isubcase in table:
                is_data = True
                case = table[isubcase]
                #print(case)
                is_real = case.is_real()
                if case.nonlinear_factor is not None:
                    times = case._times
                    is_static = False
                else:
                    is_static = True
                    times = np.zeros(1, dtype='int32')
                #print('times = ', times)
                break
                #return is_data, is_static, is_real, times
        #print('isubcase =', isubcase)
        return is_data, is_static, is_real, times

    def _get_stress_times(self, model, isubcase):
        table_types = self._get_stress_table_types()
        is_real = True
        is_data = False
        is_static = False
        times = None
        for table_type in table_types:
            if not hasattr(model, table_type):
                # print('no table_type=%s' % table_type)
                continue
            table = getattr(model, table_type)
            if isubcase in table:
                is_data = True
                case = table[isubcase]
                is_real = case.is_real()
                if case.nonlinear_factor is not None:
                    times = case._times
                    is_static = False
                else:
                    is_static = True
                    times = np.zeros(1, dtype='int32')
                break
                #return is_data, is_static, is_real, times
        return is_data, is_static, is_real, times

    def _get_stress_table_types(self):
        """
        Gets the list of Nastran stress objects that the GUI supports
        """
        table_types = [
            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress',
            'celas2_stress',
            'celas3_stress',
            'celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain',
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'crod_stress',
            'conrod_stress',
            'ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'crod_strain',
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress
            'cbar_stress',
            # OES - isotropic CBAR strain
            'cbar_strain',
            # OES - isotropic CBEAM stress
            'cbeam_stress',
            # OES - isotropic CBEAM strain
            'cbeam_strain',

            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress',
            'cquad4_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain',
            'cquad4_strain',

            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress',
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain',
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress
            'cshear_stress',
            # OES - CSHEAR strain
            'cshear_strain',
            # OES - CEALS1 224, CELAS3 225
            'nonlinear_spring_stress',
            # OES - GAPNL 86
            'nonlinear_cgap_stress',
            # OES - CBUSH 226
            'nolinear_cbush_stress',
        ]

        table_types += [
            # OES - CTRIAX6
            'ctriax_stress',
            'ctriax_strain',

            'cbush_stress',
            'cbush_strain',
            'cbush1d_stress_strain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinear_rod_stress',
            'nonlinear_rod_strain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinear_plate_stress',
            'nonlinear_plate_strain',
            #'hyperelastic_plate_stress',
            'hyperelastic_cquad4_strain',

            # OES - composite CTRIA3/CQUAD4 stress
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',

            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',

            # OGS1 - grid point stresses
            'grid_point_stresses',        # tCode=26
            'grid_point_volume_stresses',  # tCode=27
        ]
        return table_types

    def _get_nastran_header(self, case, dt, itime):
        #if case is None:
            #return None
        try:
            code_name = case.data_code['name']
        except KeyError:
            return 'Static'

        if isinstance(dt, float):
            header = ' %s = %.4E' % (code_name, dt)
        else:
            header = ' %s = %i' % (code_name, dt)

        # cases:
        #   1. lsftsfqs
        #   2. loadIDs, eigrs
        #   3. lsdvmns, eigrs
        #   ???
        if hasattr(case, 'mode_cycle'):
            header += '; freq = %g Hz' % case.mode_cycles[itime]
        elif hasattr(case, 'cycles'):
            header += '; freq = %g Hz' % case.cycles[itime]
        elif hasattr(case, 'eigis'):
            eigi = case.eigis[itime]
            cycle = np.abs(eigi) / (2. * np.pi)
            header += '; freq = %g Hz' % cycle
        elif hasattr(case, 'eigns'):  #  eign is not eigr; it's more like eigi
            eigi2 = case.eigrs[itime] #  but |eigi| = sqrt(|eign|)
            cycle = np.sqrt(np.abs(eigr2)) / (2. * np.pi)
            header += '; freq = %g Hz' % cycle
        elif hasattr(case, 'dt'):
            time = case._times[itime]
            header += 'time = %g sec' % time
            pass
        elif hasattr(case, 'lftsfqs') or hasattr(case, 'lsdvmns') or hasattr(case, 'loadIDs'):
            pass
            #raise RuntimeError(header)
        else:
            msg = 'unhandled case; header=%r\n%s' % (header, str(case))
            print(msg)
            #raise RuntimeError(msg)

        return header.strip('; ')

    def _fill_op2_time_centroidal_strain_energy(self, cases, model,
                                                key, icase, itime,
                                                form_dict, header_dict, is_static):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        case = None
        subcase_id = key[2]
        strain_energies = [
            (model.cquad4_strain_energy, 'CQUAD4', True),
            (model.cquad8_strain_energy, 'CQUAD8', True),
            (model.cquadr_strain_energy, 'CQUADR', True),
            (model.cquadx_strain_energy, 'CQUADX', True),

            (model.ctria3_strain_energy, 'CTRIA3', True),
            (model.ctria6_strain_energy, 'CTRIA6', True),
            (model.ctriar_strain_energy, 'CTRIAR', True),
            (model.ctriax_strain_energy, 'CTRIAX', True),
            (model.ctriax6_strain_energy, 'CTRIAX6', True),

            (model.ctetra_strain_energy, 'CTETRA', True),
            (model.cpenta_strain_energy, 'CPENTA', True),
            (model.chexa_strain_energy, 'CHEXA', True),
            (model.cpyram_strain_energy, 'CPYRAM', True),

            (model.crod_strain_energy, 'CROD', True),
            (model.ctube_strain_energy, 'CTUBE', True),
            (model.conrod_strain_energy, 'CONROD', True),

            (model.cbar_strain_energy, 'CBAR', True),
            (model.cbeam_strain_energy, 'CBEAM', True),

            (model.cgap_strain_energy, 'CGAP', True),
            (model.celas1_strain_energy, 'CELAS1', True),
            (model.celas2_strain_energy, 'CELAS2', True),
            (model.celas3_strain_energy, 'CELAS3', True),
            (model.celas4_strain_energy, 'CELAS4', True),
            (model.cdum8_strain_energy, 'CDUM8', False),
            (model.cbush_strain_energy, 'CBUSH', True),
            #(model.chexa8fd_strain_energy, '', False),
            (model.cbend_strain_energy, 'CBEND', False),
            (model.dmig_strain_energy, 'DMIG', False),
            (model.genel_strain_energy, 'GENEL', False),
            (model.cshear_strain_energy, 'CSHEAR', True),
        ]
        has_strain_energy = [key in res[0] for res in strain_energies]
        if not any(has_strain_energy):
            return icase
        itrue = has_strain_energy.index(True)
        ese0 = strain_energies[itrue][0]
        #times = ese0._times

        #fmt = '%g'
        #header = ''
        #form0 = ('Element Strain Energy', None, [])

        #op2.strain_energy[1]
            #type=StrainEnergyObject ntimes=3 nelements=16
            #energy, percent, density
            #modes = [1, 2, 3]

        nelements = self.nElements

        eids = self.element_ids
        ese = np.full(nelements, np.nan, dtype='float32')
        percent = np.full(nelements, np.nan, dtype='float32')
        strain_energy_density = np.full(nelements, np.nan, dtype='float32')
        for i, is_true in enumerate(has_strain_energy):
            if not is_true:
                continue
            resdict, name, flag = strain_energies[i]

            #print('key =', key)
            case = resdict[key]

            if case.is_complex():
                continue
            itotal = np.where(case.element[itime, :] == 100000000)[0][0]
            #print('itotal = ', itotal)

            eidsi2 = case.element[itime, :itotal]
            i = np.searchsorted(eids, eidsi2)
            if len(i) != len(np.unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                #print('eids = %s\n' % str(list(eids)))
                #print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)
            ese[i] = case.data[itime, :itotal, 0]
            percent[i] = case.data[itime, :itotal, 1]
            strain_energy_density[i] = case.data[itime, :itotal, 2]


        #ese

        # helicopter.dat
        #CBEAM : 10
        #CQUAD4 : 11388
        #CROD : 544
        #CTRIA3 : 151
        # nelements = 12093

        #try:
            #header = self._get_nastran_header(case, dt, itime)
            #header_dict[(key, itime)] = header
        #except AttributeError:
            #pass
        if 1:
            ese_res = GuiResult(subcase_id, header='Strain Energy',
                                title='Strain Energy', data_format='%.3e',
                                location='centroid', scalar=ese)
            percent_res = GuiResult(subcase_id, header='Percent of Total',
                                title='Percent of Total', data_format='%.3f',
                                location='centroid', scalar=percent)
            sed_res = GuiResult(subcase_id, header='Strain Energy Density',
                                title='Strain Energy Density', data_format='%.3e',
                                location='centroid', scalar=strain_energy_density)

            cases[icase] = (ese_res, (subcase_id, 'Strain Energy'))
            cases[icase + 1] = (percent_res, (subcase_id, 'Percent'))
            cases[icase + 2] = (sed_res, (subcase_id, 'Strain Energy Density'))

            form_dict[(key, itime)].append(('Strain Energy', icase, []))
            form_dict[(key, itime)].append(('Percent', icase + 1, []))
            form_dict[(key, itime)].append(('Strain Energy Density', icase + 1, []))
            icase += 3

        return icase

    #icase = self._fill_op2_time_centroidal_force(
        #cases, model, subcase_id, icase, itime, form_dict,
        #is_static)

    def _create_op2_time_centroidal_force_arrays(self, model, nelements, key, itime, header_dict):
        """
        creates the following force outputs:
         - fx, fy, fz, mx, my, mz
         - thermal_load
        """
        fx = np.zeros(nelements, dtype='float32') # axial
        fy = np.zeros(nelements, dtype='float32') # shear_y
        fz = np.zeros(nelements, dtype='float32') # shear_z

        rx = np.zeros(nelements, dtype='float32') # torque
        ry = np.zeros(nelements, dtype='float32') # bending_y
        rz = np.zeros(nelements, dtype='float32') # bending_z

        is_element_on = np.zeros(nelements, dtype='float32') # torque
        fmt = '%g'
        header = ''
        form0 = ('Force', None, [])

        case = None
        found_force = False
        for res_type in (model.conrod_force, model.crod_force, model.ctube_force):
            if key in res_type:
                found_force = True
                case = res_type[key]
                if case.is_complex():
                    continue
                data = case.data
                if case.nonlinear_factor is None:
                    ntimes = data.shape[:1]
                    eids = case.element
                    dt = case._times[itime]
                    header = self._get_nastran_header(case, dt, itime)
                    header_dict[(key, itime)] = header
                    #eids_to_find = intersect1d(self.element_ids, eids)
                    i = np.searchsorted(self.element_ids, eids)
                    assert np.array_equal(self.element_ids[i], eids)
                    fxi = data[itime, :, 0]
                    rxi = data[itime, :, 1]
                    if fxi.size != i.size:
                        msg = 'fx.size=%s i.size=%s fx=%s eids_to_find=%s' % (
                            fxi.size, i.size, fxi, eids)
                        raise RuntimeError(msg)
                    fx[i] = fxi
                    rx[i] = rxi
                    is_element_on[i] = 1.
                else:
                    continue

        if key in model.cbar_force:
            found_force = True
            ## CBAR-34
            case = model.cbar_force[key]
            if case.is_real():
                eids = case.element
                i = np.searchsorted(self.element_ids, eids)
                is_element_on[i] = 1.

                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header

                #[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2,
                # shear1, shear2, axial, torque]
                #fx[i] = case.data[:, :, 6]
                #fy[i] = case.data[:, :, 4]
                #fz[i] = case.data[:, :, 5]

                if i.size == 1:
                    rxi = case.data[itime, :, 7].max()
                    ryi = np.vstack([case.data[itime, :, 0], case.data[itime, :, 2]]).max()
                    rzi = np.vstack([case.data[itime, :, 1], case.data[itime, :, 3]]).max()
                else:
                    rxi = case.data[itime, :, 7]#.max(axis=0)
                    ryi = np.vstack([case.data[itime, :, 0], case.data[itime, :, 2]]).max(axis=0)
                    rzi = np.vstack([case.data[itime, :, 1], case.data[itime, :, 3]]).max(axis=0)
                    rzv = rzi

                    # rza = array([case.data[itime, :, 1], case.data[itime, :, 3]])#.max(axis=0)
                    # rzh = hstack([case.data[itime, :, 1], case.data[itime, :, 3]])#.max(axis=0)
                    # print(rzv.shape, rzv.shape, rzv.shape)
                assert rxi.size == i.size, 'rx.size=%s i.size=%s rx=%s' % (rxi.size, i.size, rxi)
                assert ryi.size == i.size, 'ry.size=%s i.size=%s ry=%s' % (ryi.size, i.size, ryi)
                assert rzi.size == i.size, 'rz.size=%s i.size=%s rz=%s' % (rzi.size, i.size, rzi)

                rx[i] = rxi
                ry[i] = ryi
                rz[i] = rzi

        if key in model.cbar_force_10nodes:
            found_force = True
            ## CBAR-100
            case = model.cbar_force_10nodes[key]
            eids = case.element
            ueids = np.unique(eids)

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            j = np.searchsorted(self.element_ids, ueids)
            is_element_on[j] = 1.
            di = j[1:-1] - j[0:-2]
            if di.max() != 2:
                #print('di =', np.unique(di))
                # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
                ii = 0
                eid_old = eids[0]
                fxi = defaultdict(list)
                fyi = defaultdict(list)
                fzi = defaultdict(list)
                rxi = defaultdict(list)
                ryi = defaultdict(list)
                rzi = defaultdict(list)
                for ii, eid in enumerate(eids):
                    fxi[eid].append(case.data[:, ii, 5])
                    fyi[eid].append(case.data[:, ii, 3])
                    fzi[eid].append(case.data[:, ii, 4])

                    rxi[eid].append(case.data[:, ii, 6])
                    ryi[eid].append(case.data[:, ii, 1])
                    rzi[eid].append(case.data[:, ii, 2])
                    #if eidi == eid_old:
                    #    fx[ii] = array([case.data[:, j, 5], case.data[:, j, 5]]).max(axis=0)
                    #else:
                for ii, eidi in zip(j, eids[j]):
                    fx[ii] = max(fxi[eidi])
                    fy[ii] = max(fyi[eidi])
                    fz[ii] = max(fyi[eidi])
                    rx[ii] = max(rxi[eidi])
                    ry[ii] = max(ryi[eidi])
                    rz[ii] = max(rzi[eidi])
            else:
                # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
                neids = len(np.unique(eids)) * 2
                if len(eids) != len(np.unique(eids)) * 2:
                    msg = 'CBAR-100 Error: len(eids)=%s neids=%s' % (len(eids), neids)
                    raise RuntimeError(msg)
                fx[i] = np.array(
                    [case.data[itime, ::-1, 5],
                     case.data[itime, 1::-1, 5]]).max(axis=0)
                fy[i] = np.array(
                    [case.data[itime, ::-1, 3],
                     case.data[itime, 1::-1, 3]]).max(axis=0)
                fz[i] = np.array(
                    [case.data[itime, ::-1, 4],
                     case.data[itime, 1::-1, 4]]).max(axis=0)
                rx[i] = np.array(
                    [case.data[itime, ::-1, 6],
                     case.data[itime, 1::-1, 6]]).max(axis=0)
                ry[i] = np.array(
                    [case.data[itime, ::-1, 1],
                     case.data[itime, 1::-1, 1]]).max(axis=0)
                rz[i] = np.array(
                    [case.data[itime, ::-1, 2],
                     case.data[itime, 1::-1, 2]]).max(axis=0)
        return found_force, fx, fy, fz, rx, ry, rz, is_element_on

    def _fill_op2_time_centroidal_force(self, cases, model,
                                        key, icase, itime,
                                        form_dict, header_dict, is_static):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        nelements = self.nElements
        out = self._create_op2_time_centroidal_force_arrays(model, nelements, key, itime, header_dict)
        found_force, fx, fy, fz, rx, ry, rz, is_element_on = out

        #new_cases = True
        subcase_id = key[2]
        if found_force:
            fmt = '%.4f'
            # header = self._get_nastran_header(case, dt, itime)

            #num_on = nelements
            num_off = 0
            if itime == 0 and is_element_on.min() == 0.0:
                ioff = np.where(is_element_on == 0)[0]
                num_off = len(ioff)
                print('force_eids_off = %s; n=%s' % (self.element_ids[ioff], num_off))
                self.log_error('force_eids_off = %s; n=%s' % (self.element_ids[ioff], num_off))
                force_on_res = GuiResult(subcase_id, header='Force - IsElementOn',
                                         title='Force\nIsElementOn',
                                         location='centroid', scalar=is_element_on)
                cases[icase] = (force_on_res, (subcase_id, 'Force\nIsElementOn'))
                form_dict[(key, itime)].append(('Force - IsElementOn', icase, []))
                #num_on -= num_off
                icase += 1

            if fx.min() != fx.max() or rx.min() != rx.max() and not num_off == nelements:
                fx_res = GuiResult(subcase_id, header='Axial', title='Axial',
                                   location='centroid', scalar=fx)
                fy_res = GuiResult(subcase_id, header='ShearY', title='ShearY',
                                   location='centroid', scalar=fy)
                fz_res = GuiResult(subcase_id, header='ShearZ', title='ShearZ',
                                   location='centroid', scalar=fz)
                mx_res = GuiResult(subcase_id, header='Torsion', title='Torsion',
                                   location='centroid', scalar=rx)
                my_res = GuiResult(subcase_id, header='BendingY', title='BendingY',
                                   location='centroid', scalar=ry)
                mz_res = GuiResult(subcase_id, header='BendingZ', title='BendingZ',
                                   location='centroid', scalar=rz)
                cases[icase] = (fx_res, (subcase_id, 'Axial'))
                cases[icase + 1] = (fy_res, (subcase_id, 'ShearY'))
                cases[icase + 2] = (fz_res, (subcase_id, 'ShearZ'))
                cases[icase + 3] = (mx_res, (subcase_id, 'Torsion'))
                cases[icase + 4] = (my_res, (subcase_id, 'BendingY'))
                cases[icase + 5] = (mz_res, (subcase_id, 'BendingZ'))

                form_dict[(key, itime)].append(('Axial', icase, []))
                form_dict[(key, itime)].append(('ShearY', icase + 1, []))
                form_dict[(key, itime)].append(('ShearZ', icase + 2, []))
                form_dict[(key, itime)].append(('Torque', icase + 3, []))
                form_dict[(key, itime)].append(('BendingY', icase + 4, []))
                form_dict[(key, itime)].append(('BendingZ', icase + 5, []))
                icase += 6

                is_axial = np.zeros(self.nElements, dtype='int8')
                is_shear_y = np.zeros(self.nElements, dtype='int8')
                is_shear_z = np.zeros(self.nElements, dtype='int8')
                is_torsion = np.zeros(self.nElements, dtype='int8')
                is_bending_y = np.zeros(self.nElements, dtype='int8')
                is_bending_z = np.zeros(self.nElements, dtype='int8')
                is_axial[np.where(np.abs(fx) > 0.0)[0]] = 1
                is_shear_y[np.where(np.abs(fy) > 0.0)[0]] = 1
                is_shear_z[np.where(np.abs(fz) > 0.0)[0]] = 1
                is_torsion[np.where(np.abs(rx) > 0.0)[0]] = 1
                is_bending_y[np.where(np.abs(ry) > 0.0)[0]] = 1
                is_bending_z[np.where(np.abs(rz) > 0.0)[0]] = 1
                #is_bending[where(abs(rx) > 0.0)[0]] = 1

                is_fx_res = GuiResult(subcase_id, header='IsAxial', title='IsAxial',
                                      location='centroid', scalar=is_axial, data_format=fmt)
                is_fy_res = GuiResult(subcase_id, header='IsShearY', title='IsShearY',
                                      location='centroid', scalar=is_shear_y, data_format=fmt)
                is_fz_res = GuiResult(subcase_id, header='IsShearZ', title='IsShearZ',
                                      location='centroid', scalar=is_shear_z, data_format=fmt)
                is_mx_res = GuiResult(subcase_id, header='IsTorsion', title='IsTorsion',
                                      location='centroid', scalar=is_torsion, data_format=fmt)
                is_my_res = GuiResult(subcase_id, header='IsBendingY', title='IsBendingY',
                                      location='centroid', scalar=is_bending_y, data_format=fmt)
                is_mz_res = GuiResult(subcase_id, header='IsBendingZ', title='IsBendingZ',
                                      location='centroid', scalar=is_bending_z, data_format=fmt)

                cases[icase] = (is_fx_res, (subcase_id, 'IsAxial'))
                cases[icase + 1] = (is_fy_res, (subcase_id, 'IsShearY'))
                cases[icase + 2] = (is_fz_res, (subcase_id, 'IsShearZ'))
                cases[icase + 3] = (is_mx_res, (subcase_id, 'IsTorsion'))
                cases[icase + 4] = (is_my_res, (subcase_id, 'IsBendingY'))
                cases[icase + 5] = (is_mz_res, (subcase_id, 'IsBendingZ'))

                form_dict[(key, itime)].append(('IsAxial', icase, []))
                form_dict[(key, itime)].append(('IsShearY', icase + 1, []))
                form_dict[(key, itime)].append(('IsShearZ', icase + 2, []))
                form_dict[(key, itime)].append(('IsTorsion', icase + 3, []))
                form_dict[(key, itime)].append(('IsBendingY', icase + 4, []))
                form_dict[(key, itime)].append(('IsBendingZ', icase + 5, []))
                icase += 6
        return icase

    def _fill_op2_time_centroidal_stress(self, cases, model, key, icase, itime,
                                         form_dict, header_dict, is_static,
                                         is_stress=True):
        """
        Creates the time accurate stress objects for the pyNastranGUI
        """
        #new_cases = True
        case = None
        #assert isinstance(subcase_id, int), type(subcase_id)
        assert isinstance(icase, int), icase
        #assert isinstance(itime, int), type(itime)
        assert is_stress in [True, False], is_stress
        eids = self.element_ids
        assert len(eids) > 0, eids
        nelements = self.nElements
        dt = None

        is_element_on = np.zeros(nelements, dtype='int8')  # is the element supported
        oxx = np.zeros(nelements, dtype='float32')
        oyy = np.zeros(nelements, dtype='float32')
        ozz = np.zeros(nelements, dtype='float32')

        txy = np.zeros(nelements, dtype='float32')
        tyz = np.zeros(nelements, dtype='float32')
        txz = np.zeros(nelements, dtype='float32')

        max_principal = np.zeros(nelements, dtype='float32')  # max
        mid_principal = np.zeros(nelements, dtype='float32')  # mid
        min_principal = np.zeros(nelements, dtype='float32')  # min
        ovm = np.zeros(nelements, dtype='float32')

        vm_word = None
        if is_stress:
            rods = [model.crod_stress, model.conrod_stress, model.ctube_stress,]
        else:
            rods = [model.crod_strain, model.conrod_strain, model.ctube_strain,]

        for result in rods:
            if key not in result:
                continue

            case = result[key]
            if case.is_complex():
                continue
            eidsi = case.element
            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                print('eids = %s\n' % str(list(eids)))
                print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)

            is_element_on[i] = 1
            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            # data=[1, nnodes, 4] where 4=[axial, SMa, torsion, SMt]
            oxx[i] = case.data[itime, :, 0]
            txy[i] = case.data[itime, :, 2]
            ovm[i] = np.sqrt(oxx[i]**2 + 3*txy[i]**2) # plane stress
            # max_principal[i] = sqrt(oxx[i]**2 + txy[i]**2)
            # min_principal[i] = max_principal[i] - 2 * txy[i]
            # simplification of:
            #   eig(A) = [oxx, txy]
            #            [txy, 0.0]
            # per Equation 7: http://www.soest.hawaii.edu/martel/Courses/GG303/Eigenvectors.pdf
            max_principal[i] = (oxx[i] + np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
            min_principal[i] = (oxx[i] - np.sqrt(oxx[i]**2 + 4 * txy[i]**2)) / 2.
        del rods


        if is_stress:
            bars = model.cbar_stress
        else:
            bars = model.cbar_strain

        if key in bars:
            case = bars[key]
            if case.is_complex():
                pass
            else:
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                #s1a = case.data[itime, :, 0]
                #s2a = case.data[itime, :, 1]
                #s3a = case.data[itime, :, 2]
                #s4a = case.data[itime, :, 3]

                axial = case.data[itime, :, 4]
                smaxa = case.data[itime, :, 5]
                smina = case.data[itime, :, 6]
                #MSt = case.data[itime, :, 7]

                #s1b = case.data[itime, :, 8]
                #s2b = case.data[itime, :, 9]
                #s3b = case.data[itime, :, 10]
                #s4b = case.data[itime, :, 11]

                smaxb = case.data[itime, :, 12]
                sminb = case.data[itime, :, 13]
                #MSc   = case.data[itime, :, 14]

                eidsi = case.element # [:, 0]

                i = np.searchsorted(eids, eidsi)
                if len(i) != len(np.unique(i)):
                    print('ibar = %s' % i)
                    print('eids = %s' % eids)
                    msg = 'ibar=%s is not unique' % str(i)
                    raise RuntimeError(msg)

                is_element_on[i] = 1.
                oxx[i] = axial

                ## TODO :not sure if this block is general for multiple CBAR elements
                samax = np.amax([smaxa, smaxb], axis=0)
                samin = np.amin([smaxa, smaxb], axis=0)
                assert len(samax) == len(i), len(samax)
                assert len(samin) == len(i)
                savm = np.amax(np.abs(
                    [smina, sminb,
                     smaxa, smaxb, axial]), axis=0)

                max_principal[i] = samax
                min_principal[i] = samin
                ovm[i] = savm
                del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
        del bars


        if is_stress:
            bars2 = model.cbar_stress_10nodes
        else:
            bars2 = model.cbar_strain_10nodes

        if key in bars2:
            case = bars2[key]
            if case.is_complex():
                pass
            else:
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                #  0    1    2    3    4     5     6     7     8
                # [sd, sxc, sxd, sxe, sxf, axial, smax, smin, MS]

                eidsi = case.element # [:, 0]
                ueidsi = np.unique(eidsi)
                istart = np.searchsorted(eidsi, ueidsi)
                iend = np.hstack(istart[1:], [len(eidsi)])
                axial = case.data[itime, :, 5]

                nbars = len(eidsi) // 10
                assert nbars * 10 == len(eidsi), 'nbars=%s neids=%s' % (nbars, len(eidsi))
                axial = case.data[itime, :, 5].reshape(nbars, 10).min(axis=1)
                smax = case.data[itime, :, 6].reshape(nbars, 10).max(axis=1)
                smin = case.data[itime, :, 7].reshape(nbars, 10).min(axis=1)


                i = np.searchsorted(eids, eidsi)
                if len(i) != len(np.unique(i)):
                    print('ibar = %s' % i)
                    print('eids = %s' % eids)
                    msg = 'ibar=%s is not unique' % str(i)
                    raise RuntimeError(msg)

                is_element_on[i] = 1.
                oxx[i] = axial

                ## TODO :not sure if this block is general for multiple CBAR elements
                svm = np.amax(np.abs([smin, smin]), axis=0)

                max_principal[i] = smax
                min_principal[i] = smin
                ovm[i] = svm
                #del axial, smaxa, smina, smaxb, sminb, eidsi, i, samax, samin, savm
        del bars2


        if is_stress:
            beams = model.cbeam_stress
        else:
            beams = model.cbeam_strain

        if key in beams:
            case = beams[key]
            if case.is_complex():
                pass
            else:
                eidsi = case.element_node[:, 0]
                ueids = np.unique(eidsi)
                #neids = len(ueids)

                # sxc, sxd, sxe, sxf
                # smax, smin, MSt, MSc
                dt = case._times[itime]
                header = self._get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                sxc = case.data[itime, :, 0]
                sxd = case.data[itime, :, 1]
                sxe = case.data[itime, :, 2]
                sxf = case.data[itime, :, 3]
                smax = case.data[itime, :, 4]
                smin = case.data[itime, :, 5]

                imin = np.searchsorted(eidsi, ueids)
                imax = np.searchsorted(eidsi, ueids, side='right')
                #sxxi = smax[imin:imax]
                for eid, imini, imaxi in zip(ueids, imin, imax):
                    oxxi = 0.
                    smaxi = 0.
                    smini = 0.
                    eid2 = self.eid_map[eid]
                    is_element_on[eid2] = 1.
                    oxxi = max(
                        sxc[imini:imaxi].max(),
                        sxd[imini:imaxi].max(),
                        sxe[imini:imaxi].max(),
                        sxf[imini:imaxi].max(),
                    )
                    smaxi = smax[imini:imaxi].max()
                    smini = smin[imini:imaxi].min()
                    ovmi = max(np.abs(smaxi), np.abs(smini))
                    oxxi = oxx[eid2]
                    max_principal[eid2] = smaxi
                    min_principal[eid2] = smini
                    ovm[eid2] = ovmi
                del eidsi, ueids, sxc, sxd, sxe, sxf, smax, smin, oxxi, smaxi, smini, ovmi
        del beams


        if is_stress:
            plates = [
                model.ctria3_stress, model.cquad4_stress,
                model.ctria6_stress, model.cquad8_stress,
                model.ctriar_stress, model.cquadr_stress,
            ]
        else:
            plates = [
                model.ctria3_strain, model.cquad4_strain,
                model.ctria6_strain, model.cquad8_strain,
                model.ctriar_strain, model.cquadr_strain,
            ]

        for result in plates:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            nlayers_per_element = nnodes_per_element * 2  # *2 for every other layer
            eidsi = case.element_node[::nlayers_per_element, 0]  # ::2 is for layer skipping

            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('iplate = %s' % i)
                print('eids = %s' % eids)
                print('eidsiA = %s' % case.element_node[:, 0])
                print('eidsiB = %s' % eidsi)
                msg = 'iplate=%s is not unique' % str(i)
                raise RuntimeError(msg)
            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            is_element_on[i] = 1.
            ntotal = case.data.shape[1]  # (ndt, ntotal, nresults)
            if nlayers_per_element == 1:
                j = None
            else:
                j = np.arange(ntotal)[::nlayers_per_element]

            #self.data[self.itime, self.itotal, :] = [fd, oxx, oyy,
            #                                         txy, angle,
            #                                         majorP, minorP, ovm]
            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            oxxi = case.data[itime, j, 1]
            oyyi = case.data[itime, j, 2]
            txyi = case.data[itime, j, 3]
            o1i = case.data[itime, j, 5]
            o3i = case.data[itime, j, 6]
            ovmi = case.data[itime, j, 7]

            for inode in range(1, nlayers_per_element):
                #print('%s - ilayer = %s' % (case.element_name, inode))
                oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 1]]), axis=0)
                oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 5]]), axis=0)
                o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 6]]), axis=0)
                ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 7]]), axis=0)
                assert len(oxxi) == len(j)

            oxx[i] = oxxi
            oyy[i] = oyyi
            txy[i] = txyi
            max_principal[i] = o1i
            min_principal[i] = o3i
            ovm[i] = ovmi


        if is_stress:
            cplates = [
                ('CTRIA3', model.ctria3_composite_stress),
                ('CQUAD4', model.cquad4_composite_stress),
                ('CTRIA6', model.ctria6_composite_stress),
                ('CQUAD8', model.cquad8_composite_stress),
                #model.ctriar_composite_stress,
                #model.cquadr_composite_stress,
            ]
        else:
            cplates = [
                ('CTRIA3', model.ctria3_composite_strain),
                ('CQUAD4', model.cquad4_composite_strain),
                ('CTRIA6', model.ctria6_composite_strain),
                ('CQUAD8', model.cquad8_composite_strain),
                #model.ctriar_composite_strain,
                #model.cquadr_composite_strain,
            ]

        for cell_type, result in cplates:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            eidsi = case.element_layer[:, 0]
            layers = case.element_layer[:, 1]
            ntotal = case.data.shape[1]

            #[o11, o22, t12, t1z, t2z, angle, major, minor, max_shear]
            oxxs = case.data[itime, :, 0]
            oyys = case.data[itime, :, 1]
            txys = case.data[itime, :, 2]
            txzs = case.data[itime, :, 3]
            tyzs = case.data[itime, :, 4]
            # angle
            omaxs = case.data[itime, :, 6]
            omins = case.data[itime, :, 7]
            ovms = case.data[itime, :, 8]

            j = 0
            for eid in np.unique(eidsi):
                ieid = np.where(eidsi == eid)[0]
                ieid.sort()
                layersi = layers[ieid]
                eid2 = self.eid_map[eid]
                is_element_on[eid2] = 1.

                oxxi = 0.
                oyyi = 0.
                txyi = 0.
                tyzi = 0.
                txzi = 0.
                omaxi = 0.
                omini = 0.
                ovmi = 0.
                nlayers = len(layersi)
                for ilayer in range(nlayers):
                    oxxi = max(oxxs[j], oxxi)
                    oyyi = max(oyys[j], oyyi)
                    txyi = max(txys[j], txyi)
                    tyzi = max(tyzs[j], tyzi)
                    txzi = max(txzs[j], txzi)

                    omaxi = max(omaxs[j], omaxi)
                    omini = min(omins[j], omini)
                    ovmi = max(ovms[j], ovmi)
                    j += 1

                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                txy[eid2] = txyi
                tyz[eid2] = tyzi
                txz[eid2] = txzi
                max_principal[eid2] = omaxi
                min_principal[eid2] = omini
                ovm[eid2] = ovmi
            del oxxi, oyyi, txyi, tyzi, txzi, omaxi, omini, ovmi, eid2, j, layers, eidsi
        del cplates


        if is_stress:
            solids = [(model.ctetra_stress),
                      (model.cpenta_stress),
                      (model.chexa_stress),]
        else:
            solids = [(model.ctetra_strain),
                      (model.cpenta_strain),
                      (model.chexa_strain),]

        for result in solids:
            if key not in result:
                continue
            case = result[key]
            if case.is_complex():
                continue

            if case.is_von_mises():
                vm_word = 'vonMises'
            else:
                vm_word = 'maxShear'

            nnodes_per_element = case.nnodes
            eidsi = case.element_cid[:, 0]
            ntotal = len(eidsi)  * nnodes_per_element

            i = np.searchsorted(eids, eidsi)
            if len(i) != len(np.unique(i)):
                print('isolid = %s' % str(i))
                print('eids = %s' % eids)
                print('eidsi = %s' % eidsi)
                assert len(i) == len(np.unique(i)), 'isolid=%s is not unique' % str(i)

            is_element_on[i] = 1
            #self.data[self.itime, self.itotal, :] = [oxx, oyy, ozz,
            #                                         txy, tyz, txz,
            #                                         o1, o2, o3, ovm]

            if nnodes_per_element == 1:
                j = None
            else:
                j = np.arange(ntotal)[::nnodes_per_element]
                ueidsi = np.unique(eidsi)
                assert len(j) == len(ueidsi), 'j=%s ueidsi=%s' % (j, ueidsi)

            dt = case._times[itime]
            header = self._get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            oxxi = case.data[itime, j, 0]
            oyyi = case.data[itime, j, 1]
            ozzi = case.data[itime, j, 2]
            txyi = case.data[itime, j, 3]
            tyzi = case.data[itime, j, 4]
            txzi = case.data[itime, j, 5]
            o1i = case.data[itime, j, 6]
            o2i = case.data[itime, j, 7]
            o3i = case.data[itime, j, 8]
            ovmi = case.data[itime, j, 9]

            for inode in range(1, nnodes_per_element):
                oxxi = np.amax(np.vstack([oxxi, case.data[itime, j + inode, 0]]), axis=0)
                oyyi = np.amax(np.vstack([oyyi, case.data[itime, j + inode, 1]]), axis=0)
                ozzi = np.amax(np.vstack([ozzi, case.data[itime, j + inode, 2]]), axis=0)
                txyi = np.amax(np.vstack([txyi, case.data[itime, j + inode, 3]]), axis=0)
                tyzi = np.amax(np.vstack([tyzi, case.data[itime, j + inode, 4]]), axis=0)
                txzi = np.amax(np.vstack([txzi, case.data[itime, j + inode, 2]]), axis=0)

                o1i = np.amax(np.vstack([o1i, case.data[itime, j + inode, 6]]), axis=0)
                o2i = np.amax(np.vstack([o2i, case.data[itime, j + inode, 7]]), axis=0)
                o3i = np.amin(np.vstack([o3i, case.data[itime, j + inode, 8]]), axis=0)
                ovmi = np.amax(np.vstack([ovmi, case.data[itime, j + inode, 9]]), axis=0)
                assert len(oxxi) == len(j)

            oxx[i] = oxxi
            oyy[i] = oyyi
            ozz[i] = ozzi
            txy[i] = txyi
            tyz[i] = tyzi
            txz[i] = txzi
            max_principal[i] = o1i
            mid_principal[i] = o2i
            min_principal[i] = o3i
            ovm[i] = ovmi
        del solids


        if is_stress:
            word = 'Stress'
            fmt = '%.3f'
        else:
            word = 'Strain'
            fmt = '%.4e'

        # a form is the table of output...
        # Subcase 1         <--- formi  - form_isubcase
        #    Time 1
        #        Stress     <--- form0  - the root level
        #            oxx    <--- formis - form_itime_stress
        #            oyy
        #            ozz

        if dt is None:
            return icase

        header = ''
        if not is_static:
            #print('is_static = %s' % is_static)
            if case is None:
                formis = None
                return icase
            header = self._get_nastran_header(case, dt, itime)
            #form_time[0] = header

        form0 = (word, None, [])
        formis = form0[2]
        # subcase_id, icase, resultType, vector_size, location, dataFormat
        subcase_id = key[2]
        if is_stress and itime == 0:
            if is_element_on.min() == 0:  # if all elements aren't on
                ioff = np.where(is_element_on == 0)[0]
                print('stress_eids_off = %s' % np.array(self.element_ids[ioff]))
                self.log_error('stress_eids_off = %s' % self.element_ids[ioff])
                stress_res = GuiResult(
                    subcase_id, header='Stress - isElementOn', title='Stress\nisElementOn',
                    location='centroid', scalar=oxx, data_format=fmt)
                cases[icase] = (stress_res, (subcase_id, 'Stress - isElementOn'))
                form_dict[(key, itime)].append(('Stress - IsElementOn', icase, []))
                icase += 1

        if oxx.min() != oxx.max():
            oxx_res = GuiResult(subcase_id, header=word + 'XX', title=word + 'XX',
                                location='centroid', scalar=oxx, data_format=fmt)
            cases[icase] = (oxx_res, (subcase_id, word + 'XX'))
            form_dict[(key, itime)].append((word + 'XX', icase, []))
            icase += 1
        if oyy.min() != oyy.max():
            oyy_res = GuiResult(subcase_id, header=word + 'YY', title=word + 'YY',
                                location='centroid', scalar=oyy, data_format=fmt)
            cases[icase] = (oyy_res, (subcase_id, word + 'YY'))
            form_dict[(key, itime)].append((word + 'YY', icase, []))
            icase += 1
        if ozz.min() != ozz.max():
            ozz_res = GuiResult(subcase_id, header=word + 'ZZ', title=word + 'ZZ',
                                location='centroid', scalar=ozz, data_format=fmt)
            cases[icase] = (ozz_res, (subcase_id, word + 'ZZ'))
            form_dict[(key, itime)].append((word + 'ZZ', icase, []))
            icase += 1
        if txy.min() != txy.max():
            oxy_res = GuiResult(subcase_id, header=word + 'XY', title=word + 'XY',
                                location='centroid', scalar=txy, data_format=fmt)
            cases[icase] = (oxy_res, (subcase_id, word + 'XY'))
            form_dict[(key, itime)].append((word + 'XY', icase, []))
            icase += 1
        if tyz.min() != tyz.max():
            oyz_res = GuiResult(subcase_id, header=word + 'YZ', title=word + 'YZ',
                                location='centroid', scalar=tyz, data_format=fmt)
            cases[icase] = (oyz_res, (subcase_id, word + 'YZ'))
            form_dict[(key, itime)].append((word + 'YZ', icase, []))
            icase += 1
        if txz.min() != txz.max():
            oxz_res = GuiResult(subcase_id, header=word + 'XZ', title=word + 'XZ',
                                location='centroid', scalar=txz, data_format=fmt)
            cases[icase] = (oxz_res, (subcase_id, word + 'XZ'))
            form_dict[(key, itime)].append((word + 'XZ', icase, []))
            icase += 1
        if max_principal.min() != max_principal.max():
            maxp_res = GuiResult(subcase_id, header='MaxPrincipal', title='MaxPrincipal',
                                 location='centroid', scalar=max_principal, data_format=fmt)
            cases[icase] = (maxp_res, (subcase_id, 'MaxPrincipal'))
            form_dict[(key, itime)].append(('Max Principal', icase, []))
            icase += 1
        if mid_principal.min() != mid_principal.max():
            midp_res = GuiResult(subcase_id, header='MidPrincipal', title='MidPrincipal',
                                 location='centroid', scalar=mid_principal, data_format=fmt)
            cases[icase] = (midp_res, (subcase_id, 'MidPrincipal'))
            form_dict[(key, itime)].append(('Mid Principal', icase, []))
            icase += 1
        if min_principal.min() != min_principal.max():
            minp_res = GuiResult(subcase_id, header='MinPrincipal', title='MinPrincipal',
                                 location='centroid', scalar=min_principal, data_format=fmt)
            cases[icase] = (minp_res, (subcase_id, 'MinPrincipal'))
            form_dict[(key, itime)].append(('Min Principal', icase, []))
            icase += 1
        if vm_word is not None:
            if not is_stress:
                max_min = max(ovm.max(), np.abs(ovm.min()))
                if max_min > 100:
                    raise RuntimeError('vm strain = %s' % ovm)

            ovm_res = GuiResult(subcase_id, header=vm_word, title=vm_word,
                                location='centroid', scalar=ovm, data_format=fmt)
            cases[icase] = (ovm_res, (subcase_id, vm_word))
            form_dict[(key, itime)].append((vm_word, icase, []))
            icase += 1

        #, case, header, form0
        return icase


    def _get_nastran_key_order(self, model):
        displacement_like = [
            model.displacements,
            model.velocities,
            model.accelerations,
            model.eigenvectors,
            model.spc_forces,
            model.mpc_forces,

            # untested
            model.load_vectors,
            model.applied_loads,
            model.force_vectors,
            #[model.grid_point_forces, 'GridPointForces'],  # TODO: this is buggy...
        ]
        temperature_like = [
            model.temperatures,
        ]
        stress = [
            model.crod_stress, model.conrod_stress, model.ctube_stress,
            model.cbar_stress,
            model.cbeam_stress,

            model.ctria3_stress, model.cquad4_stress,
            model.ctria6_stress, model.cquad8_stress,
            model.ctriar_stress, model.cquadr_stress,

            model.ctria3_composite_stress, model.cquad4_composite_stress,
            model.ctria6_composite_stress, model.cquad8_composite_stress,

            model.ctetra_stress, model.cpenta_stress, model.chexa_stress,
        ]
        strain = [
            model.crod_strain, model.conrod_strain, model.ctube_strain,
            model.cbar_strain,
            model.cbeam_strain,

            model.ctria3_strain, model.cquad4_strain,
            model.ctria6_strain, model.cquad8_strain,
            model.ctriar_strain, model.cquadr_strain,

            model.ctria3_composite_strain, model.cquad4_composite_strain,
            model.ctria6_composite_strain, model.cquad8_composite_strain,

            model.ctetra_strain, model.cpenta_strain, model.chexa_strain,
        ]
        strain_energy = [
            #model.strain_energy,
            model.cquad4_strain_energy, model.cquad8_strain_energy,
            model.cquadr_strain_energy, model.cquadx_strain_energy,
            model.ctria3_strain_energy, model.ctria6_strain_energy,
            model.ctriar_strain_energy, model.ctriax_strain_energy,
            model.ctriax6_strain_energy,
            model.ctetra_strain_energy, model.cpenta_strain_energy,
            model.chexa_strain_energy,
            model.crod_strain_energy, model.ctube_strain_energy,
            model.conrod_strain_energy,
            model.cbar_strain_energy, model.cbeam_strain_energy,
            model.cgap_strain_energy, model.cbush_strain_energy,
            model.celas1_strain_energy, model.celas2_strain_energy,
            model.celas3_strain_energy, model.celas4_strain_energy,
            model.cdum8_strain_energy, model.dmig_strain_energy,
            model.cbend_strain_energy,
            model.genel_strain_energy, model.cshear_strain_energy,
        ]

        result_groups = [
            displacement_like, temperature_like, stress, strain, strain_energy,
        ]

        nids = self.node_ids
        eids = self.element_ids
        keys_order = []
        # model = OP2()

        # subcase_ids = model.subcase_key.keys()

        #self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.analysis_code, self.label]
        subcase_ids = list(model.iSubcaseNameMap.keys())
        subcase_ids.sort()
        #print('subcase_ids =', subcase_ids)


        # isubcase, analysis_code, sort_method, count, subtitle
        #(1, 2, 1, 0, 'SUPERELEMENT 0') : result1

        #subcase_ids = subcase_ids
        #print('subcase_idsB =' % subcase_ids)
        for isubcase in sorted(subcase_ids):
            if isubcase == 0:
                # beam_modes
                self.log.error('*****isubcase=0')
                continue
            # value = (analysis_codei, sort_methodi, counti, isubtitle)
            #print('subcase_key =', model.subcase_key)
            keys = model.subcase_key[isubcase]
            #print('keys[%s] =%s' % (isubcase, keys))
            key0 = keys[0]

            # this while loop lets us make sure we pull the analysis codes in the expected order
            # TODO: doesn't pull count in the right order
            # TODO: doesn't pull subtitle in right order
            keys2 = deepcopy(keys)
            while keys2:
                key = keys2[-1]
                #print('while keys ->', key)
                (analysis_code, sort_method, count, subtitle) = key
                #assert isubcase == isubcasei, 'isubcase=%s isubcasei=%s' % (isubcase, isubcasei)
                assert analysis_code < 12, analysis_code
                for ianalysis_code in range(12):
                    keyi = (ianalysis_code, sort_method, count, subtitle)
                    if keyi in keys2:
                        #print(keyi)
                        keyi2 = (isubcase, ianalysis_code, sort_method, count, subtitle)
                        #print(keyi2)
                        keys_order.append(keyi2)
                        keys2.remove(keyi)
                #keys_order += keys
        return keys_order
