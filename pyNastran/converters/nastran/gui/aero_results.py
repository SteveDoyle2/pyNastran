from __future__ import annotations
from typing import TYPE_CHECKING
from itertools import count
import numpy as np

from pyNastran.utils import PathLike
from pyNastran.bdf.bdf import BDF

from pyNastran.op2.op2 import OP2
from pyNastran.op2.result_objects.stress_object import _get_nastran_header
from pyNastran.gui.gui_objects.gui_result import GuiResult
from pyNastran.gui.gui_objects.force_results import ForceResults2
from pyNastran.gui.gui_objects.displacement_results import DisplacementResults2
from .result_objects.displacement_reduced import DisplacementReduced
from .result_objects.case2d import Case2D

if TYPE_CHECKING:
    from cpylog import SimpleLogger
    from .nastran_io import NastranIO


def is_early_return_aero(self: NastranIO, model: OP2) -> bool:
    """identify an aero model"""
    if not self.aero_is_quad_mesh:
        return False

    early_return_aero = False
    # for aero identification
    nnode = len(self.node_ids)
    nelement = len(self.element_ids)

    for key, case in model.displacements.items():
        ngrid = len(case.node_gridtype)
        if ngrid > nnode:
            # I hate nastran aero; this isn't robust.
            # it puts [nodes, aero_nodes, aero_elements] in the result.
            # This is broken by adding a single grid in the model.
            early_return_aero = True
            break
    for key, case in model.eigenvectors.items():
        ngrid = len(case.node_gridtype)
        if ngrid > nnode:
            # I hate nastran aero; this isn't robust.
            # it puts [nodes, aero_nodes, aero_elements] in the result.
            # this is broken by adding a single grid in the model.
            early_return_aero = True
            break

    trim_results = model.op2_results.trim
    if trim_results.aero_pressure:
        # model.log.error(f'fem: nnode={nnode} nelement={nelement}')
        for key, case in trim_results.aero_pressure.items():
            # case.cp
            # case.pressure
            # case.nodes
            ncp = len(case.cp)
            if ncp == nelement:
                early_return_aero = True
                break
    elif trim_results.aero_force:
        for key, case in trim_results.aero_force.items():
            # aero_force[8]:
            #   nodes:        n=3108
            #   force:        (3108, 6)
            #   force_label:  (3108,)
            # print(case.get_stats())
            nforce = len(case.force)
            if nforce == nelement:
                early_return_aero = True
                break
    return early_return_aero

def load_nastran_results_aero(results_filename: PathLike,
                              model_aero: BDF,
                              results_model: OP2,
                              xyz_cid0: np.ndarray,
                              cases, form, icase: int):
    """
    create results for aero models
     - displacements or real/complex eigenvectors
     - no spc_forces or load_vectors for aero

    TODO: improve case tagging (e.g., mode=3; freq=10 Hz) in the lower-left corner
    """
    aero_nids = np.array(list(model_aero.nodes), dtype='int32')
    aero_eids = np.array(list(model_aero.elements), dtype='int32')

    subcase_keys = list(results_model.displacements)
    trim_results = results_model.op2_results.trim
    for key in results_model.eigenvectors:
        if key not in subcase_keys:
            subcase_keys.append(key)
    for key in trim_results.aero_pressure:
        if key not in subcase_keys:
            subcase_keys.append(key)
    for key in trim_results.aero_force:
        if key not in subcase_keys:
            subcase_keys.append(key)

    subcase_id_to_subcase_word = {}
    mydicts = [
        results_model.displacements,
        results_model.eigenvectors,
        trim_results.aero_pressure,
        trim_results.aero_force,
    ]
    for key in subcase_keys:
        for mydict in mydicts:
            if key not in mydict:
                continue
            subcase_id = key[0]
            if subcase_id in subcase_id_to_subcase_word:
                continue
            case = mydict[key]
            subtitle = case.subtitle
            label = case.label
            subcase_word = f'Subcase {subcase_id}'
            if subtitle:
                subcase_word += f'; subtitle={subtitle}'
            if subtitle:
                subcase_word += f'; label={label}'
            subcase_id_to_subcase_word[subcase_id] = subcase_word

    results_filename = str(results_filename).strip(r'.\\')
    results_form = []
    mmax = xyz_cid0.max(axis=0)
    mmin = xyz_cid0.min(axis=0)
    dim_max = (mmax - mmin).max()
    for key in subcase_keys:
        # key = (8, 1, 1, 0, 0, '', '')
        # (3, 2, 1, 0, 0, '', '')
        # adding aero eigenvector (real)
        # (3, 9, 1, 0, 0, '', '')
        # adding aero eigenvector (complex)

        subcase_form = []
        subcase_id = key[0]
        try:
            subcase_word = subcase_id_to_subcase_word[subcase_id]
        except KeyError:
            print(f'key={key}')
            print(subcase_id_to_subcase_word)
            raise
        if key in results_model.displacements:
            icase = _aero_deflection(
                results_model.displacements,
                results_model.log,
                aero_nids,
                aero_eids,
                xyz_cid0, dim_max,
                cases, subcase_form, icase, key,
                resname='Deflection')
        if key in results_model.eigenvectors:
            icase = _aero_deflection(
                results_model.eigenvectors,
                results_model.log,
                aero_nids,
                aero_eids,
                xyz_cid0, dim_max,
                cases, subcase_form, icase, key,
                resname='Eigenvector')

        if key in trim_results.aero_pressure:
            # print(f'key = {key}')
            aero_pressure = trim_results.aero_pressure[key]
            cp = aero_pressure.cp
            cp_res = GuiResult(
                subcase_id, header='Aero Cp', title='Aero Cp',
                location='centroid', scalar=cp)
            # (8, 1, 1, 0, 0, '', '')
            # print(aero_pressure.get_stats())
            cases[icase] = (cp_res, (0, 'Aero Cp'))
            # Subcase 1: Aero Cp

            # self.title = title
            # self.subtitle = subtitle
            # self.label = label
            subcase_form.append(('Aero Cp', icase, []))
            icase += 1

        if key in trim_results.aero_force:
            aero_force = trim_results.aero_force[key]
            # aero_force[8]:
            # fem: nnode=4522 nelement=3108
            #   nodes:        n=3108
            #   force:        (3108, 6)
            #   force_label:  (3108,)
            fz_res = GuiResult(
                subcase_id, header='Aero Force - Fz', title='Aero Force - Fz',
                location='centroid', scalar=aero_force.force[:, 2])
            my_res = GuiResult(
                subcase_id, header='Aero Force - My', title='Aero Force - My',
                location='centroid', scalar=aero_force.force[:, 4])

            # ---- elemental_forces_to_nodal_forces ----
            nodal_force_dict = {}
            nodal_moment_dict = {}
            nnodes_dict = {}
            for nid in aero_nids:
                nodal_force_dict[nid] = np.zeros(3)
                nodal_moment_dict[nid] = np.zeros(3)
                nnodes_dict[nid] = 0
            for eid, forcei in zip(aero_eids, aero_force.force):
                elem = model_aero.elements[eid]
                element_normal = elem.Normal()
                for nid in elem.nodes:
                    nodal_force_dict[nid] += forcei[2] * element_normal  # Fz
                    nodal_moment_dict[nid] += forcei[4] * element_normal  # My
                    nnodes_dict[nid] += 1

            naero_node = len(aero_nids)
            nodal_forces = np.zeros((naero_node, 3), dtype='float64')
            nodal_moments = np.zeros((naero_node, 3), dtype='float64')
            for i, (nid, nnodei) in zip(count(), nnodes_dict.items()):
                force = nodal_force_dict[nid]
                moment = nodal_moment_dict[nid]
                nnodei = nnodes_dict[nid]
                nodal_forces[i] = force / nnodei
                nodal_moments[i] = moment / nnodei
            #-----------------------------------------------

            force_case = Case2D(aero_nids, nodal_forces)
            moment_case = Case2D(aero_nids, nodal_moments)
            methods_txyz_rxyz_force = ['Fx', 'Fy', 'Fz']
            methods_txyz_rxyz_moment = ['Mx', 'My', 'Mz']
            index_to_base_title_annotation_force = {
                0: {'title': 'F_', 'corner': 'F_'},
            }
            index_to_base_title_annotation_moment = {
                0: {'title': 'M_', 'corner': 'M_'},
            }
            force_res = ForceResults2(
                subcase_id,
                aero_nids, xyz_cid0,
                force_case, aero_force.title,
                index_to_base_title_annotation=index_to_base_title_annotation_force,
                t123_offset=0, methods_txyz_rxyz=methods_txyz_rxyz_force,
                dim_max=1.0, data_format='%g',
                is_variable_data_format=False,
                nlabels=None, labelsize=None, ncolors=None, colormap='',
                set_max_min=False, uname='NastranGeometry-ForceResults2')
            moment_res = ForceResults2(
                subcase_id,
                aero_nids, xyz_cid0,
                moment_case, aero_force.title,
                index_to_base_title_annotation=index_to_base_title_annotation_moment,
                t123_offset=0, methods_txyz_rxyz=methods_txyz_rxyz_moment,
                dim_max=1.0, data_format='%g',
                is_variable_data_format=False,
                nlabels=None, labelsize=None, ncolors=None, colormap='',
                set_max_min=False, uname='NastranGeometry-ForceResults2')
            cases[icase] = (fz_res, (0, 'Aero Force - Fz'))
            cases[icase+1] = (my_res, (0, 'Aero Force - My'))
            cases[icase+2] = (force_res, (0, 'Aero Force'))
            cases[icase+3] = (moment_res, (0, 'Aero Moment'))
            subcase_form.append(('Aero Force - Fz', icase, []))
            subcase_form.append(('Aero Force - My', icase+1, []))
            subcase_form.append(('Aero Force', icase+2, []))
            subcase_form.append(('Aero Moment', icase+3, []))
            icase += 4

        if len(subcase_form):
            results_form.append((subcase_word, None, subcase_form))

    if len(results_form):
        form.append((results_filename, None, results_form))
        # form_results = (basename + '-Results', None, form_optimization)


def _aero_deflection(displacements_dict: dict,
                     log: SimpleLogger,
                     aero_nids: np.ndarray,
                     aero_eids: np.ndarray,
                     xyz_cid0: np.ndarray,
                     dim_max: float,
                     cases, form, icase: int,
                     key: tuple,
                     resname='Deflection',
                     ) -> int:
    """
    aero deflection results are appended to the OUGV1 displacement table as:
      - [structure_nodes, aero_nodes, aero_elements]
    """
    assert resname in ['Deflection', 'Eigenvector']
    node_ids = aero_nids
    naero_nodes = len(aero_nids)
    naero_eids = len(aero_eids)

    nextra_nodes = naero_nodes + naero_eids
    log.debug(f'naero_nodes = {naero_nodes}')
    log.debug(f'naero_eids = {naero_eids}')
    log.debug(f'ntotal = {nextra_nodes}')

    # --------------------------------------------------------------------------
    log.debug('--------------------------------------------------------------')
    disp_case = displacements_dict[key]

    subcase_id = key[0]
    node_gridtype = disp_case.node_gridtype
    all_nids = node_gridtype[:, 0]
    max_nid = all_nids.max()
    nids_aero = all_nids[-naero_nodes:]
    log.debug(f'nids_aero = {nids_aero}')
    log.debug(f'max_nid = {max_nid}')

    ntimes = disp_case.data.shape[0]
    node_gridtype = disp_case.node_gridtype[-nextra_nodes:-naero_eids, :].copy()
    node_gridtype[:, 0] = node_ids
    nodal_disp = disp_case.data[:, -nextra_nodes:-naero_eids, :]
    element_disp = disp_case.data[:, -naero_eids:, :]
    log.debug(f'aero_data.shape = {str(nodal_disp.shape)}')
    log.debug(f'element_disp.shape = {str(element_disp.shape)}')
    assert nodal_disp.shape[1] == naero_nodes, nodal_disp.shape
    assert element_disp.shape[1] == naero_eids, element_disp.shape

    # ---------------------------------------------------------------------

    disp_case_aero = DisplacementReduced(disp_case, nodal_disp, node_gridtype)
    deflection_res = DisplacementResults2(
        subcase_id, node_ids, xyz_cid0, disp_case_aero,
        title=resname,
        t123_offset=0,
        dim_max=dim_max,
        data_format='%g', nlabels=None, labelsize=None,
        ncolors=None, colormap='', set_max_min=False,
        # uname=resname,
    )
    # cases[icase] = (deflection_res, (0, 'Aero Deflection'))
    # form.append(('Aero Deflection', icase, []))

    # 1: statics
    # 2: real modes
    # 9: complex modes
    is_real = (disp_case.analysis_code in {1, 2})
    is_complex = (disp_case.analysis_code == 9)
    assert is_real or is_complex, (disp_case.data_code, disp_case.get_stats())

    for itime in range(ntimes):
        # mode = 2; freq = 75.9575 Hz
        dt = disp_case._times[itime]
        header = _get_nastran_header(disp_case, dt, itime)
        # header_dict[(key, itime)] = header
        # keys_map[key] = KeyMap(disp_case.subtitle, disp_case.label,
        #                        disp_case.superelement_adaptivity_index,
        #                        disp_case.pval_step)

        # headers2.append(header)
        # cases[icase] = (deflection_res, (itime, title1))  # do I keep this???
        # formii = (title1, icase, [])
        # form_dict[(key, itime)].append(formii)
        wordi = f'Aero {resname}'
        if resname == 'Eigenvector':
            wordi += f'; {header}'
        cases[icase] = (deflection_res, (itime, wordi))
        form.append((wordi, icase, []))
        icase += 1
    return icase
