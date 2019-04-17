"""
Interface for converting OP2 results to the GUI format
"""
# pylint: disable=C1801, C0103
from __future__ import print_function
from copy import deepcopy
from collections import defaultdict

import numpy as np
from numpy.linalg import norm  # type: ignore

from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.converters.nastran.geometry_helper import NastranGuiAttributes
from pyNastran.converters.nastran.displacements import (
    DisplacementResults, ForceTableResults) #, TransientElementResults
from pyNastran.op2.result_objects.stress_object import (
    _get_nastran_header,
    get_rod_stress_strain,
    get_bar_stress_strain, get_bar100_stress_strain, get_beam_stress_strain,
    get_plate_stress_strain, get_solid_stress_strain)
from pyNastran.gui.gui_objects.gui_result import GridPointForceResult


class NastranGuiResults(NastranGuiAttributes):
    """
    Defines OP2 specific methods NastranIO
    """
    def __init__(self):
        super(NastranGuiResults, self).__init__()


    def _fill_grid_point_forces(self, cases, model, key, icase,
                                form_dict, header_dict, keys_map):
        if key not in model.grid_point_forces:
            return icase
        grid_point_forces = model.grid_point_forces[key]
        case = grid_point_forces
        if not case.is_real:
            #raise RuntimeError(grid_point_forces.is_real)
            return icase

        subcase_id = key[0]
        title = 'Grid Point Forces'
        header = 'Grid Point Forces'
        nastran_res = GridPointForceResult(subcase_id, header, title, grid_point_forces)

        itime = 0

        cases[icase] = (nastran_res, (itime, 'Grid Point Forces'))
        formii = ('Grid Point Forces', icase, [])
        form_dict[(key, itime)].append(formii)

        dt = case._times[itime]
        header = _get_nastran_header(case, dt, itime)
        header_dict[(key, itime)] = header
        keys_map[key] = (case.subtitle, case.label,
                         case.superelement_adaptivity_index, case.pval_step)

        icase += 1
        return icase

    def _fill_op2_oug_oqg(self, cases, model, key, icase,
                          form_dict, header_dict, keys_map):
        """
        loads nodal results bector results (e.g., dispalcements/temperatures)
        """
        icase = self._fill_nastran_displacements(cases, model, key, icase,
                                                 form_dict, header_dict, keys_map)
        icase = self._fill_nastran_temperatures(cases, model, key, icase,
                                                form_dict, header_dict, keys_map)
        return icase

    def _fill_nastran_displacements(self, cases, model, key, icase,
                                    form_dict, header_dict, keys_map):
        """
        loads the nodal dispalcements/velocity/acceleration/eigenvector/spc/mpc forces
        """
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

        for (result, name, deflects) in displacement_like:
            if key not in result:
                continue
            for t123_offset in [0, 3]:
                #if t123_offset == 3:
                    #continue
                try:
                    icase = self._fill_nastran_ith_displacement(
                        result, name, deflects, t123_offset,
                        cases, model, key, icase,
                        form_dict, header_dict, keys_map)
                except ValueError:
                    if not t123_offset == 3:
                        raise
                    self.log.error('skipping %s result; t123_offset=%s; type=%s' % (
                        name, t123_offset, result[key].__class__.__name__))
        return icase

    def _fill_nastran_ith_displacement(self, result, name, deflects, t123_offset,
                                       cases, model, key, icase,
                                       form_dict, header_dict, keys_map):
        """helper for ``_fill_nastran_displacements`` to unindent the code a bit"""
        nnodes = self.nnodes
        nids = self.node_ids
        if t123_offset == 0:
            title1 = name + ' T_XYZ'
        else:
            assert t123_offset == 3, t123_offset
            title1 = name + ' R_XYZ'
        #title2 = name + ' R_XYZ'

        case = result[key]
        subcase_idi = case.isubcase
        if not hasattr(case, 'data'):
            print('str(%s) has no data...' % case.__class.__name__)
            return icase
        #if not case.is_real:
            #print('complex results not supported...')
            #continue
        # transient
        if case.nonlinear_factor is not None:
            #code_name = case.data_code['name']
            unused_has_cycle = hasattr(case, 'mode_cycle')
        else:
            unused_has_cycle = False
            unused_code_name = None
        if not case.is_sort1:
            self.log.warning('Skipping because SORT2\n' + str(case))
            return icase

        t123, tnorm, ntimes = _get_t123_tnorm(case, nids, nnodes,
                                              t123_offset=t123_offset)

        titles = []
        scales = []
        headers = []
        #if deflects:
        if deflects:
            nastran_res = DisplacementResults(subcase_idi, titles, headers,
                                              self.xyz_cid0, t123, tnorm,
                                              scales,
                                              uname='NastranResult')

            dmax = []
            for itime in range(ntimes):
                dt = case._times[itime]

                if name == 'Displacement':
                    # (6673, )
                    normiii = np.linalg.norm(t123[itime, :, :], axis=1)
                    #print(normiii.shape)
                    #print('Displacement; itime=%s time=%s tnorm=%s' % (
                        #itime, dt, normiii.max()))
                    dmax.append(normiii.max())
                # mode = 2; freq = 75.9575 Hz
                header = _get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)

                tnorm_abs_max = tnorm.max()
                #if tnorm_abs_max == 0.0:
                    #scale = self.displacement_scale_factor
                #else:
                    #scale = self.displacement_scale_factor / tnorm_abs_max

                scale = self.gui.settings.dim_max
                if tnorm_abs_max > 0.0:
                    scale = self.gui.settings.dim_max / tnorm_abs_max * 0.25
                scales.append(scale)
                titles.append(title1)
                headers.append(header)
                cases[icase] = (nastran_res, (itime, title1))  # do I keep this???
                formii = (title1, icase, [])
                form_dict[(key, itime)].append(formii)
                icase += 1

            if name == 'Displacement':
                # Displacement; itime=361 time=3.61 tnorm=1.46723
                #print('dmax = ', max(dmax))
                pass
            nastran_res.save_defaults()
        else:
            nastran_res = ForceTableResults(subcase_idi, titles, headers,
                                            t123, tnorm,
                                            scales, #deflects=deflects,
                                            uname='NastranResult')
            for itime in range(ntimes):
                dt = case._times[itime]
                header = _get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)

                tnorm_abs_max = tnorm.max()
                scale = 1.
                scales.append(scale)
                titles.append(title1)
                headers.append(header)
                cases[icase] = (nastran_res, (itime, title1))  # do I keep this???
                formii = (title1, icase, [])
                form_dict[(key, itime)].append(formii)
                icase += 1
            nastran_res.save_defaults()
        return icase

    def _fill_nastran_temperatures(self, cases, model, key, icase,
                                   form_dict, header_dict, keys_map):
        """loads the nodal temperatures"""
        nnodes = self.nnodes
        #nids = self.node_ids
        temperature_like = [
            (model.temperatures, 'Temperature'),
        ]
        for (result, name) in temperature_like:
            if key not in result:
                continue
            case = result[key]
            subcase_idi = case.isubcase
            if not hasattr(case, 'data'):
                continue

            if not case.is_sort1:
                self.log.warning('Skipping because SORT2\n' + str(case))
                continue
            assert case.is_sort1, case.is_sort1

            ntimes = case.ntimes
            for itime in range(ntimes):
                dt = case._times[itime]
                header = _get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)

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

    def _fill_op2_force(self, cases, model, key, icase, itime,
                        form_dict, header_dict, keys_map):
        """creates the force plots"""
        #assert isinstance(key, int), key
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        icase = self._fill_op2_time_centroidal_force(
            cases, model, key, icase, itime,
            form_dict, header_dict, keys_map)
        return icase

    def _fill_op2_stress(self, cases, model, key, icase, itime,
                         form_dict, header_dict, keys_map, is_stress=True):
        """creates the stress plots"""
        assert isinstance(icase, int), icase
        assert isinstance(form_dict, dict), form_dict
        icase = self._fill_op2_time_centroidal_stress(
            cases, model, key, icase, itime, form_dict, header_dict, keys_map,
            is_stress=is_stress)
        return icase

    def _fill_op2_strain(self, cases, model, key, icase, itime,
                         form_dict, header_dict, keys_map):
        """creates the strain plots"""
        return self._fill_op2_stress(cases, model, key, icase, itime,
                                     form_dict, header_dict, keys_map,
                                     is_stress=False)

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
                is_real = case.is_real
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
            'grid_point_surface_stresses',        # tCode=26
            'grid_point_volume_stresses',  # tCode=27
        ]
        return table_types

    def _fill_op2_time_gpstress(self, cases, model,
                                key, icase, itime,
                                form_dict, header_dict, keys_map):
        """
        Creates the time accurate grid point stress objects for the pyNastranGUI
        """
        #print(key, icase, itime)
        if key in model.grid_point_stresses_volume_direct:
            case = model.grid_point_stresses_volume_direct[key]

            #print(''.join(case.get_stats()))
            if case.is_complex:
                return icase

            dt = case._times[itime]
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header

            # volume direct
            #['ox', 'oy', 'oz', 'txy', 'tyz', 'txz', 'pressure', 'ovm']
            nids = self.node_ids
            nnodes = len(nids)
            ox = np.full(nnodes, np.nan, dtype='float32')
            oy = np.full(nnodes, np.nan, dtype='float32')
            oz = np.full(nnodes, np.nan, dtype='float32')
            txy = np.full(nnodes, np.nan, dtype='float32')
            tyz = np.full(nnodes, np.nan, dtype='float32')
            txz = np.full(nnodes, np.nan, dtype='float32')
            ovm = np.full(nnodes, np.nan, dtype='float32')

            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)
            subcase_id = key[0]

            nids2 = case.node
            i = np.searchsorted(nids, nids2)
            if len(i) != len(np.unique(i)):
                msg = 'irod=%s is not unique\n' % str(i)
                #print('eids = %s\n' % str(list(eids)))
                #print('eidsi = %s\n' % str(list(eidsi)))
                raise RuntimeError(msg)
            ox[i] = case.data[itime, :, 0]
            oy[i] = case.data[itime, :, 1]
            oz[i] = case.data[itime, :, 2]
            txy[i] = case.data[itime, :, 3]
            tyz[i] = case.data[itime, :, 4]
            txz[i] = case.data[itime, :, 5]
            ovm[i] = case.data[itime, :, 7]

            headers = ['oxx', 'oyy', 'ozz', 'txy', 'tyz', 'txz', 'ovm']
            form = [('Volume Direct', None, [])]
            formi = form[0][2]
            form_dict[(key, itime)] = form

            for header, resi in zip(headers, (ox, oy, oz, txy, tyz, txz, ovm)):
                ese_res = GuiResult(subcase_id, header=header,
                                    title=header, data_format='%.3e',
                                    location='node', scalar=resi)
                cases[icase] = (ese_res, (subcase_id, header))
                formi.append((header, icase, []))
                icase += 1

        return icase

    def _fill_op2_time_centroidal_strain_energy(self, cases, model,
                                                key, icase, itime,
                                                form_dict, unused_header_dict, keys_map):
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
            (model.conm2_strain_energy, 'CONM2', False),
        ]
        has_strain_energy = [key in res[0] for res in strain_energies]
        if not any(has_strain_energy):
            return icase
        itrue = has_strain_energy.index(True)
        unused_ese0 = strain_energies[itrue][0]
        #times = ese0._times

        #fmt = '%g'
        #header = ''
        #form0 = ('Element Strain Energy', None, [])

        #op2.strain_energy[1]
            #type=StrainEnergyObject ntimes=3 nelements=16
            #energy, percent, density
            #modes = [1, 2, 3]

        nelements = self.nelements

        eids = self.element_ids
        ese = np.full(nelements, np.nan, dtype='float32')
        percent = np.full(nelements, np.nan, dtype='float32')
        strain_energy_density = np.full(nelements, np.nan, dtype='float32')
        for i, is_true in enumerate(has_strain_energy):
            if not is_true:
                continue
            resdict, unused_name, unused_flag = strain_energies[i]

            #print('key =', key)
            case = resdict[key]
            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)

            if case.is_complex:
                continue
            itotal = np.where(case.element[itime, :] == 100000000)[0][0]
            #print('itotal = ', itotal)

            eidsi2 = case.element[itime, :itotal]
            i = np.searchsorted(eids, eidsi2)
            if len(i) != len(np.unique(i)):
                msg = 'i%s=%s is not unique' % (case.element_name, str(i))
                #print('eids = %s\n' % str(list(eids)))
                #print('eidsi = %s\n' % str(list(eidsi)))
                model.log.warning(msg)
                continue
                #raise RuntimeError(msg)

            # verifies the try-except is what we think it is (missing elements)
            esei = case.data[itime, :itotal, 0]
            try:
                ese[i] = esei
                percent[i] = case.data[itime, :itotal, 1]
                strain_energy_density[i] = case.data[itime, :itotal, 2]
            except IndexError:
                model.log.warning('error reading Strain Energy')
                continue

        #ese

        # helicopter.dat
        #CBEAM : 10
        #CQUAD4 : 11388
        #CROD : 544
        #CTRIA3 : 151
        # nelements = 12093

        #try:
            #header = _get_nastran_header(case, dt, itime)
            #header_dict[(key, itime)] = header
        #except AttributeError:
            #pass

        if np.any(np.isfinite(ese)):
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

    def _create_op2_time_centroidal_force_arrays(self, model, nelements, key, itime,
                                                 header_dict, keys_map):
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
        unused_fmt = '%g'
        header = ''
        unused_form0 = ('Force', None, [])

        case = None
        found_force = False
        for res_type in (model.conrod_force, model.crod_force, model.ctube_force):
            if key in res_type:
                found_force = True
                case = res_type[key]
                if case.is_complex:
                    continue
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)
                data = case.data
                if case.nonlinear_factor is None:
                    unused_ntimes = data.shape[:1]
                    eids = case.element
                    dt = case._times[itime]
                    header = _get_nastran_header(case, dt, itime)
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
            if case.is_real:
                eids = case.element
                i = np.searchsorted(self.element_ids, eids)
                is_element_on[i] = 1.

                dt = case._times[itime]
                header = _get_nastran_header(case, dt, itime)
                header_dict[(key, itime)] = header
                keys_map[key] = (case.subtitle, case.label,
                                 case.superelement_adaptivity_index, case.pval_step)

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
                    unused_rzv = rzi

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
            header = _get_nastran_header(case, dt, itime)
            header_dict[(key, itime)] = header
            keys_map[key] = (case.subtitle, case.label,
                             case.superelement_adaptivity_index, case.pval_step)

            j = np.searchsorted(self.element_ids, ueids)
            di = j[1:-1] - j[0:-2]
            if len(di) == 0:
                # pload1
                self.log_error('Error loading CBAR-100 forces; failed slicing element_ids')
            else:
                is_element_on[j] = 1.

                if di.max() != 2:
                    #print('di =', np.unique(di))
                    # [station, bending_moment1, bending_moment2, shear1, shear2, axial, torque]
                    ii = 0
                    unused_eid_old = eids[0]
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
                                        form_dict, header_dict, keys_map):
        """
        Creates the time accurate strain energy objects for the pyNastranGUI
        """
        nelements = self.nelements
        out = self._create_op2_time_centroidal_force_arrays(
            model, nelements, key, itime, header_dict, keys_map)
        found_force, fx, fy, fz, rx, ry, rz, is_element_on = out

        #new_cases = True
        subcase_id = key[2]
        if found_force:
            fmt = '%.4f'
            # header = _get_nastran_header(case, dt, itime)

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

                is_axial = np.zeros(self.nelements, dtype='int8')
                is_shear_y = np.zeros(self.nelements, dtype='int8')
                is_shear_z = np.zeros(self.nelements, dtype='int8')
                is_torsion = np.zeros(self.nelements, dtype='int8')
                is_bending_y = np.zeros(self.nelements, dtype='int8')
                is_bending_z = np.zeros(self.nelements, dtype='int8')
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
                                         form_dict, header_dict, keys_map,
                                         is_stress=True):
        """
        Creates the time accurate stress objects for the pyNastranGUI
        """
        #new_cases = True
        #assert isinstance(subcase_id, int), type(subcase_id)
        assert isinstance(icase, int), icase
        #assert isinstance(itime, int), type(itime)
        assert is_stress in [True, False], is_stress
        eids = self.element_ids
        assert len(eids) > 0, eids
        nelements = self.nelements

        is_element_on = np.zeros(nelements, dtype='int8')  # is the element supported
        oxx = np.full(nelements, np.nan, dtype='float32')
        oyy = np.full(nelements, np.nan, dtype='float32')
        ozz = np.full(nelements, np.nan, dtype='float32')

        txy = np.full(nelements, np.nan, dtype='float32')
        tyz = np.full(nelements, np.nan, dtype='float32')
        txz = np.full(nelements, np.nan, dtype='float32')

        max_principal = np.full(nelements, np.nan, dtype='float32')  # max
        mid_principal = np.full(nelements, np.nan, dtype='float32')  # mid
        min_principal = np.full(nelements, np.nan, dtype='float32')  # min
        #max_shear = np.full(nelements, np.nan, dtype='float32')
        ovm = np.full(nelements, np.nan, dtype='float32')

        vm_word = None
        #-------------------------------------------------------------
        #vm_word = get_spring_stress_strain(
            #model, key, is_stress, vm_word, itime,
            #oxx, txy,
            #max_principal, min_principal, ovm, is_element_on,
            #eids, header_dict, keys_map)

        #-------------------------------------------------------------
        vm_word = get_rod_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx, txy,
            max_principal, min_principal, ovm, is_element_on,
            eids, header_dict, keys_map)

        vm_word = get_bar_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx,
            max_principal, min_principal, ovm, is_element_on,
            eids, header_dict, keys_map)

        vm_word = get_bar100_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx,
            max_principal, min_principal, ovm, is_element_on,
            eids, header_dict, keys_map)

        vm_word = get_beam_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx,
            max_principal, min_principal, ovm, is_element_on,
            header_dict, keys_map, self.eid_map)
        #-------------------------------------------------------------
        vm_word = get_plate_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx, oyy, txy, max_principal, min_principal, ovm, is_element_on,
            eids, header_dict, keys_map)

        #vm_word = get_shear_stress_strain(
            #model, key, is_stress, vm_word, itime,
            #oxx, txy,
            #max_principal, min_principal, ovm, is_element_on,
            #eids, header_dict, keys_map)

        if is_stress:
            stress_obj = self.stress[key]
        else:
            stress_obj = self.strain[key]

        if len(stress_obj.composite_data_dict):
            str(stress_obj)
            vm_word = stress_obj.set_composite_stress_old(
                key, itime, oxx, oyy, txy, tyz, txz,
                max_principal, min_principal, ovm,
                is_element_on, header_dict,
            )

        vm_word = get_solid_stress_strain(
            model, key, is_stress, vm_word, itime,
            oxx, oyy, ozz, txy, tyz, txz,
            max_principal, mid_principal, min_principal, ovm, is_element_on,
            eids, header_dict, keys_map)

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

        if vm_word is None:
            #print('vm_word is None')
            return icase

        form0 = (word, None, [])
        unused_formis = form0[2]
        subcase_id = key[2]
        if is_stress and itime == 0:
            if is_element_on.min() == 0:  # if all elements aren't on
                print_empty_elements(self.model, self.element_ids, is_element_on, self.log_error)

                stress_res = GuiResult(
                    subcase_id, header='Stress - isElementOn', title='Stress\nisElementOn',
                    location='centroid', scalar=oxx, data_format=fmt)
                cases[icase] = (stress_res, (subcase_id, 'Stress - isElementOn'))
                form_dict[(key, itime)].append(('Stress - IsElementOn', icase, []))
                icase += 1

        #print('max/min', max_principal.max(), max_principal.min())
        if np.any(np.isfinite(oxx)):
            oxx_res = GuiResult(subcase_id, header=word + 'XX', title=word + 'XX',
                                location='centroid', scalar=oxx, data_format=fmt)
            cases[icase] = (oxx_res, (subcase_id, word + 'XX'))
            form_dict[(key, itime)].append((word + 'XX', icase, []))
            icase += 1

        if np.any(np.isfinite(oyy)):
            oyy_res = GuiResult(subcase_id, header=word + 'YY', title=word + 'YY',
                                location='centroid', scalar=oyy, data_format=fmt)
            cases[icase] = (oyy_res, (subcase_id, word + 'YY'))
            form_dict[(key, itime)].append((word + 'YY', icase, []))
            icase += 1

        if np.any(np.isfinite(ozz)):
            ozz_res = GuiResult(subcase_id, header=word + 'ZZ', title=word + 'ZZ',
                                location='centroid', scalar=ozz, data_format=fmt)
            cases[icase] = (ozz_res, (subcase_id, word + 'ZZ'))
            form_dict[(key, itime)].append((word + 'ZZ', icase, []))
            icase += 1

        if np.any(np.isfinite(txy)):
            oxy_res = GuiResult(subcase_id, header=word + 'XY', title=word + 'XY',
                                location='centroid', scalar=txy, data_format=fmt)
            cases[icase] = (oxy_res, (subcase_id, word + 'XY'))
            form_dict[(key, itime)].append((word + 'XY', icase, []))
            icase += 1

        if np.any(np.isfinite(tyz)):
            oyz_res = GuiResult(subcase_id, header=word + 'YZ', title=word + 'YZ',
                                location='centroid', scalar=tyz, data_format=fmt)
            cases[icase] = (oyz_res, (subcase_id, word + 'YZ'))
            form_dict[(key, itime)].append((word + 'YZ', icase, []))
            icase += 1

        if np.any(np.isfinite(txz)):
            oxz_res = GuiResult(subcase_id, header=word + 'XZ', title=word + 'XZ',
                                location='centroid', scalar=txz, data_format=fmt)
            cases[icase] = (oxz_res, (subcase_id, word + 'XZ'))
            form_dict[(key, itime)].append((word + 'XZ', icase, []))
            icase += 1

        if np.any(np.isfinite(max_principal)):
            maxp_res = GuiResult(subcase_id, header='MaxPrincipal', title='MaxPrincipal',
                                 location='centroid', scalar=max_principal, data_format=fmt)
            cases[icase] = (maxp_res, (subcase_id, 'MaxPrincipal'))
            form_dict[(key, itime)].append(('Max Principal', icase, []))
            icase += 1

        if np.any(np.isfinite(mid_principal)):
            midp_res = GuiResult(subcase_id, header='MidPrincipal', title='MidPrincipal',
                                 location='centroid', scalar=mid_principal, data_format=fmt)
            cases[icase] = (midp_res, (subcase_id, 'MidPrincipal'))
            form_dict[(key, itime)].append(('Mid Principal', icase, []))
            icase += 1

        if np.any(np.isfinite(min_principal)):
            minp_res = GuiResult(subcase_id, header='MinPrincipal', title='MinPrincipal',
                                 location='centroid', scalar=min_principal, data_format=fmt)
            cases[icase] = (minp_res, (subcase_id, 'MinPrincipal'))
            form_dict[(key, itime)].append(('Min Principal', icase, []))
            icase += 1

        if vm_word is not None:
            ovm_res = GuiResult(subcase_id, header=vm_word, title=vm_word,
                                location='centroid', scalar=ovm, data_format=fmt)
            cases[icase] = (ovm_res, (subcase_id, vm_word))
            form_dict[(key, itime)].append((vm_word, icase, []))
            icase += 1

        #, case, header, form0
        return icase

def print_empty_elements(model, element_ids, is_element_on, log_error):
    """prints the first 20 elements that aren't supportedas part of the stress results"""
    ioff = np.where(is_element_on == 0)[0]
    print('stress_eids_off = %s' % np.array(element_ids[ioff]))
    log_error('stress_eids_off = %s' % element_ids[ioff])

    i = 0
    imax = 20
    for eid in element_ids[ioff]:
        element = model.elements[eid]
        if element.type not in ['CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4']:
            print(element.rstrip())
            i += 1
            if i == imax:
                break
    print('-----------------------------------')


def _get_t123_tnorm(case, nids, nnodes, t123_offset=0):
    """
    helper method for _fill_op2_oug_oqg

    Parameters
    ----------
    case : DisplacementArray, ForceArray, etc.
        the OP2 result object???
    nids : (nnodes,) int ndarray
        the nodes in the model???
    nnodes : int
        the number of nodes in the model???
    t123_offset : int; default=0
        0 : translations / forces
        3 : rotations / moments

    """
    assert case.is_sort1, case.is_sort1

    itime0 = 0
    t1 = case.data[itime0, :, 0]
    ndata = t1.shape[0]
    if nnodes != ndata:
        #print('nnodes=%s ndata=%s' % (nnodes, ndata))
        nidsi = case.node_gridtype[:, 0]
        #assert len(nidsi) == nnodes, 'nidsi=%s nnodes=%s' % (nidsi, nnodes)
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
    t123 = case.data[:, :, t123_offset:t123_offset+3]
    ntimes = case.ntimes

    if nnodes != ndata:
        t123i = np.zeros((ntimes, nnodes, 3), dtype='float32')
        t123i[:, j, :] = t123
        t123 = t123i

        # (itime, nnodes, xyz)
        # tnorm (901, 3)
        tnorm = norm(t123, axis=2)   # I think this is wrong...
        #print('tnorm.shape ', tnorm.shape)
        assert len(tnorm) == t123.shape[0]
    else:
        # (itime, nnodes, xyz)
        # tnorm (901, 3)

        # float32s are apparently buggy in numpy if you have small numbers
        # see models/elements/loadstep_elememnts.op2
        try:
            tnorm = norm(t123, axis=1)
        except FloatingPointError:
            t123 = t123.astype(dtype='float64')
            tnorm = norm(t123, axis=1)

            #print('skipping %s' % name)
            #print(t123.max(axis=1))
            #for itime, ti in enumerate(t123):
                #print('itime=%s' % itime)
                #print(ti.tolist())
        assert len(tnorm) == t123.shape[0]

    assert t123.shape[0] == ntimes, 'shape=%s expected=(%s, %s, 3)' % (t123.shape, ntimes, nnodes)
    assert t123.shape[1] == nnodes, 'shape=%s expected=(%s, %s, 3)' % (t123.shape, ntimes, nnodes)
    return t123, tnorm, ntimes


def _get_times(model, key):
    """
    Get the times/frequencies/eigenvalues/loadsteps used on a given
    subcase
    """
    table_types = model.get_table_types()
    is_real = True
    is_data = False
    is_static = False
    times = None
    for table_type in table_types:
        if not model.has_result(table_type):
            #model.log.debug('no table_type=%s' % table_type)
            continue
        table = model.get_result(table_type)
        if len(table) == 0:
            continue
        #print(key, table, type(table))

        if key in table:
            is_data = True
            case = table[key]
            #print(case)
            is_real = case.is_real

            # you're presumably looking here because of a bug
            # are you sure the keys are the right length?
            #print("is_real=%r nonlinear_factor=%r _times=%s" % (
                #is_real, case.nonlinear_factor, case._times))
            if case.nonlinear_factor is not None:
                times = case._times
                is_static = False
            else:
                is_static = True
                times = np.zeros(1, dtype='int32')
            #print('times = ', times)
            break
            #return is_data, is_static, is_real, times
    return is_data, is_static, is_real, times
