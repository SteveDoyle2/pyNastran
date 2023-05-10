"""tests the NastranIO class"""
import os
from typing import Set

import vtk
from cpylog import SimpleLogger

import pyNastran
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.gui.nastran_io import NastranIO

from pyNastran.gui.gui_objects.gui_result import GuiResult, NormalResult, GridPointForceResult
from pyNastran.gui.gui_objects.displacements import ForceTableResults, DisplacementResults, ElementalTableResults
from pyNastran.converters.nastran.gui.results import SimpleTableResults, LayeredTableResults



class NastranGUI(NastranIO, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)
        self.build_fmts(['nastran'], stop_on_failure=True)


def _save_nastran_results(test: NastranGUI) -> vtkUnstructuredGrid:
    log = test.log
    vtk_ugrid = test.grid

    point_data = vtk_ugrid.GetPointData()
    cell_data = vtk_ugrid.GetCellData()

    used_titles = set()
    for key, case_data in test.result_cases.items():
        case, index_name = case_data
        if case.is_complex:
            log.warning(f'skipping case {str(case)} because it is complex')
            continue
        if isinstance(case, (GridPointForceResult)):
            log.warning(f'skipping case {str(case)}')
            continue

        if isinstance(case, ForceTableResults):
            _save_force_table_results(key, index_name, case, used_titles,
                                      point_data, cell_data, log)
            continue

        elif isinstance(case, DisplacementResults):
            _save_displacement_results(key, index_name, case, used_titles,
                                       point_data, cell_data, log)
            continue

        elif isinstance(case, SimpleTableResults):
            _save_simple_table_results(key, index_name, case, used_titles,
                                       point_data, cell_data, log)
            continue
        elif isinstance(case, LayeredTableResults):
            vtk_array = _save_layered_table_results(key, index_name, case, used_titles,
                                                    point_data, cell_data, log)
        else:
            if case.title == 'Normals' and case.scalar is None:
                continue

            titlei = case.title
            if case.subcase_id > 0:
                titlei = f'{case.title}_subcase={case.subcase_id:d}'

            vtk_array = numpy_to_vtk(case.scalar, deep=0, array_type=None)
            if titlei in used_titles:
                log.warning(f'skipping GuiResult {titlei} because it is already used')
                continue
            _check_title(titlei, used_titles)
            vtk_array.SetName(titlei)
        _add_array(case.location, point_data, cell_data, vtk_array)
    return vtk_ugrid

def _save_force_table_results(key: int, index_name,
                              case: ForceTableResults,
                              used_titles: set[str],
                              point_data, cell_data, log: SimpleLogger) -> None:
    if index_name[0] != 0:
        return
    title = case.titles[0]
    dxyz = case.dxyz

    if dxyz.ndim == 2:
        vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
        assert title not in used_titles, title
        titlei =  f'{title}_subcase={case.subcase_id}'
        _check_title(title, used_titles)
        vtk_array.SetName(titlei)
        _add_array(case.location, point_data, cell_data, vtk_array)
    elif dxyz.ndim == 3:
        for itime, title in enumerate(case.titles):
            header = case.headers[itime].replace(' = ', '=')
            dxyz = case.dxyz[itime, :, :]
            vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
            titlei =  f'{header}_subcase={case.subcase_id}'
            _check_title(titlei, used_titles)
            vtk_array.SetName(titlei)
            _add_array(case.location, point_data, cell_data, vtk_array)
    else:
        log.warning(f'cannot add {str(case)}')

def _save_displacement_results(key: int, index_name,
                               case: DisplacementResults,
                               used_titles: set[str],
                               point_data, cell_data, log: SimpleLogger) -> None:
    if index_name[0] != 0:
        return
    assert case.dxyz.ndim == 3, case.dxyz.shape
    for itime, title in enumerate(case.titles):
        #'Displacement T_XYZ: Static'
        #'Eigenvectors T_XYZ: mode = 1; freq = 9.31323e-10 Hz'
        header = case.headers[itime].replace(' = ', '=')

        dxyz = case.dxyz[itime, :, :]
        assert dxyz.ndim == 2, dxyz.shape
        vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
        titlei =  f'{header}_subcase={case.subcase_id}'
        _check_title(titlei, used_titles)
        vtk_array.SetName(titlei)
        _add_array(case.location, point_data, cell_data, vtk_array)

def _save_simple_table_results(key: int, index_name,
                               case: SimpleTableResults,
                               used_titles: set[str],
                               point_data, cell_data, log: SimpleLogger) -> None:
    # (1, (0, 0, 'Static'))
    # SimpleTableResults:
    #     title=None
    #     subcase_id=1
    #     data_type='<f4'
    #     is_real=True is_complex=False
    #     location='centroid'
    #     header=[]
    #     methods=['σxx', 'MS_axial', 'τxy', 'MS_torsion']
    #     data_format=['%.3f', '%.3f', '%.3f', '%.3f']
    #     uname='Rod Stress'

    #(itime, imethod, unused_header) = name
    #ntimes = self.scalars.shape[0]
    #j = ntimes * imethod + itime
    #(itime, imethod, unused_header) = name
    #scalars = self.scalars[itime, :, imethod]
    name = index_name[1]
    (itime, imethod, header) = name
    header2 = header.replace(' = ', '=')
    method = case.methods[imethod]
    titlei =  f'{case.uname}: {method}_{header2}_subcase={case.subcase_id}'
    res = case.get_result(key, name)

    # form_name = case.form_names[itime, ilayer, imethod]
    _check_title(titlei, used_titles)
    vtk_array = numpy_to_vtk(res, deep=0, array_type=None)
    vtk_array.SetName(titlei)

    del name, itime, imethod, header, header2
    _add_array(case.location, point_data, cell_data, vtk_array)
    #log.warning(f'skipping SimpleTableResults {case}')

def _save_layered_table_results(key: int, index_name,
                                case: LayeredTableResults,
                                used_titles: set[str],
                                point_data, cell_data, log: SimpleLogger) -> vtk.vtkFloatArray:
        assert case.location == 'centroid', case
        #(1, (0, 0, 'Static')) = index_name
        name = index_name[1]
        #print(case, name)
        (itime, ilayer, imethod, unused_header) = name

        #ntimes, nlayers, nmethods, nheader = case.scalars.shape
        #k = t0 + nmethods * ilayer + imethod
        #method = case.methods[imethod]
        #form_index = case.get_form_index(key, name)
        form_name = case.form_names[itime, ilayer, imethod]
        #for method in case.methods:
        titlei =  f'{form_name}_subcase={case.subcase_id}'
        res = case.get_result(key, name)

        _check_title(titlei, used_titles)
        vtk_array = numpy_to_vtk(res, deep=0, array_type=None)
        vtk_array.SetName(titlei)
        del name, itime, ilayer, imethod
        return vtk_array

def nastran_to_vtk(op2_filename: str, vtk_filename: str):
    test = NastranGUI()
    test.create_secondary_actors = False

    log = test.log
    log.level = 'error'
    log.level = 'warning'

    test.load_nastran_geometry(op2_filename)
    test.load_nastran_results(op2_filename)
    vtk_ugrid = _save_nastran_results(test)

    #root = vtk.vtkMultiBlockDataSet()
    #coords_branch = vtk.vtkMultiBlockDataSet()

    #root.SetBlock(0, vtk_ugrid)
    #root.SetBlock(1, coords_branch)

    #for coord_id, axes_actor in test.axes.items():
        #coords_branch.SetBlock(0, axes)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(vtk_filename)
    writer.SetInputData(vtk_ugrid)
    writer.Write()
    #print('done')


def _add_array(location, point_data, cell_data, vtk_array):
    if location == 'node':
        point_data.AddArray(vtk_array)
    else:
        assert location == 'centroid'
        cell_data.AddArray(vtk_array)


def _check_title(title: str, used_titles: set[str]) -> None:
    #print(title)
    assert title not in used_titles, title
    used_titles.add(title)

def main() -> None:
    PKG_PATH = pyNastran.__path__[0]
    MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

    op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
    vtk_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.vtu')
    nastran_to_vtk(op2_filename, vtk_filename)

    op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.op2')
    vtk_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.vtu')
    nastran_to_vtk(op2_filename, vtk_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
    print("done")
