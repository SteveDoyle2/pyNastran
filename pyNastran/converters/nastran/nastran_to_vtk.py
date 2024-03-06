"""tests the NastranIO class"""
import os

#from vtk import vtkPointData, vtkCellData, vtkFloatArray, vtkXMLUnstructuredGridWriter
from vtkmodules.vtkCommonDataModel import vtkPointData, vtkCellData
from vtkmodules.vtkCommonCore import vtkFloatArray
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter

from cpylog import SimpleLogger

import pyNastran
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk
from pyNastran.gui.testing_methods import FakeGUIMethods
from pyNastran.converters.nastran.gui.nastran_io import NastranIO

from pyNastran.gui.gui_objects.gui_result import GuiResult, GridPointForceResult, check_title # NormalResult,
from pyNastran.gui.gui_objects.displacements import ForceTableResults, DisplacementResults # , ElementalTableResults
from pyNastran.converters.nastran.gui.result_objects.simple_table_results import SimpleTableResults
from pyNastran.converters.nastran.gui.result_objects.layered_table_results import LayeredTableResults
from pyNastran.converters.nastran.gui.result_objects.displacement_results import DisplacementResults2
from pyNastran.converters.nastran.gui.result_objects.force_results import ForceResults2
from pyNastran.converters.nastran.gui.result_objects.solid_stress_results import SolidStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.composite_stress_results import CompositeStrainStressResults2
from pyNastran.converters.nastran.gui.result_objects.plate_stress_results import PlateStrainStressResults2


class NastranGUI(NastranIO, FakeGUIMethods):
    def __init__(self, inputs=None):
        FakeGUIMethods.__init__(self, inputs=inputs)
        NastranIO.__init__(self)
        self.build_fmts(['nastran'], stop_on_failure=True)


def save_nastran_results(gui: NastranGUI,
                         vtk_ugrid: vtkUnstructuredGrid) -> None:
    log = gui.log

    point_data = vtk_ugrid.GetPointData()
    cell_data = vtk_ugrid.GetCellData()

    used_titles = set()
    for key, case_data in gui.result_cases.items():
        icase = key
        case, index_name = case_data
        if case.is_complex:
            log.warning(f'skipping case {str(case)} because it is complex')
            continue
        if isinstance(case, (GridPointForceResult)):
            log.warning(f'skipping case {str(case)}')
            continue

        if isinstance(case, ForceTableResults):
            _save_force_table_results(case, key, index_name, used_titles,
                                      point_data, cell_data, log)
            continue

        elif isinstance(case, DisplacementResults):
            _save_displacement_results(case, key, index_name, used_titles,
                                       point_data, cell_data, log)
            continue

        elif isinstance(case, SimpleTableResults):
            _save_simple_table_results(case, key, index_name, used_titles,
                                       point_data, cell_data, log)
            continue

        if isinstance(case, LayeredTableResults):
            vtk_array = _save_layered_table_results(case,
                                                    key, index_name, used_titles,
                                                    point_data, cell_data, log)
            location = case.location

        elif isinstance(case, (DisplacementResults2, ForceResults2)):
            #vector = case.get_force_vector_result(*index_name)
            vector, *unused_junk = case.get_vector_data_dense(*index_name)
            titlei = f'icase={icase}; ' + case.get_annotation(*index_name).replace(' = ', '=')
            check_title(titlei, used_titles)
            vtk_array = numpy_to_vtk(vector, deep=0, array_type=None)
            vtk_array.SetName(titlei)
            location = case.location
        elif isinstance(case, (CompositeStrainStressResults2, PlateStrainStressResults2, SolidStrainStressResults2)):
            fringe = case.get_fringe_result(*index_name)
            titlei = f'icase={icase}; ' + case.get_annotation(*index_name).replace(' = ', '=')
            check_title(titlei, used_titles)
            vtk_array = numpy_to_vtk(fringe, deep=0, array_type=None)
            vtk_array.SetName(titlei)
            location = case.get_location(*index_name)
        else:
            if case.title == 'Normals' and case.scalar is None:
                continue
            if not isinstance(case, GuiResult):  # pragma: no cover
                log.warning(f'skipping case {str(case)} because it is unhandled')
                log.warning(str(case))
                return
            assert isinstance(case, GuiResult), case
            location = case.location
            vtk_array = case.save_vtk_result(used_titles)
        add_vtk_array(location, point_data, cell_data, vtk_array)

def _save_force_table_results(case: ForceTableResults,
                              key: int,
                              index_name: tuple[int, tuple[int, int, str]],
                              used_titles: set[str],
                              point_data: vtkPointData, cell_data: vtkCellData,
                              log: SimpleLogger) -> None:
    if index_name[0] != 0:
        return
    title = case.titles[0]
    dxyz = case.dxyz

    if dxyz.ndim == 2:
        vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
        assert title not in used_titles, title
        titlei =  f'{title}_subcase={case.subcase_id}'
        check_title(title, used_titles)
        vtk_array.SetName(titlei)
        add_vtk_array(case.location, point_data, cell_data, vtk_array)
    elif dxyz.ndim == 3:
        for itime, title in enumerate(case.titles):
            header = case.headers[itime].replace(' = ', '=')
            dxyz = case.dxyz[itime, :, :]
            vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
            titlei =  f'{header}_subcase={case.subcase_id}'
            check_title(titlei, used_titles)
            vtk_array.SetName(titlei)
            add_vtk_array(case.location, point_data, cell_data, vtk_array)
    else:
        log.warning(f'cannot add {str(case)!r}')

def _save_displacement_results(case: DisplacementResults,
                               key: int,
                               index_name: tuple[int, tuple[int, int, str]],
                               used_titles: set[str],
                               point_data: vtkPointData, cell_data: vtkCellData,
                               log: SimpleLogger) -> None:
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
        check_title(titlei, used_titles)
        vtk_array.SetName(titlei)
        add_vtk_array(case.location, point_data, cell_data, vtk_array)

def _save_simple_table_results(case: SimpleTableResults,
                               key: int,
                               index_name: tuple[int, tuple[int, int, str]],
                               used_titles: set[str],
                               point_data: vtkPointData, cell_data: vtkCellData,
                               log: SimpleLogger) -> None:
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
    fringe, vector = case.get_fringe_vector_result(key, name)
    res = vector if vector is not None else fringe

    # form_name = case.form_names[itime, ilayer, imethod]
    check_title(titlei, used_titles)
    vtk_array = numpy_to_vtk(res, deep=0, array_type=None)
    vtk_array.SetName(titlei)

    del name, itime, imethod, header, header2
    add_vtk_array(case.location, point_data, cell_data, vtk_array)
    #log.warning(f'skipping SimpleTableResults {case}')

def _save_layered_table_results(case: LayeredTableResults,
                                key: int,
                                index_name: tuple[int, tuple[int, int, str]],
                                used_titles: set[str],
                                point_data: vtkPointData, cell_data: vtkCellData,
                                log: SimpleLogger) -> vtkFloatArray:
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
    method = case.get_methods(key, name)[0]
    fringe, vector = case.get_fringe_vector_result(key, name)
    res = vector if vector is not None else fringe

    check_title(titlei, used_titles)
    vtk_array = numpy_to_vtk(res, deep=0, array_type=None)
    vtk_array.SetName(titlei)
    del name, itime, ilayer, imethod
    return vtk_array

def nastran_to_vtk(bdf_filename: str,
                   op2_filename: str,
                   vtk_filename: str) -> None:
    """kind of a hack, but it will always work assuming the GUI works"""
    gui = NastranGUI()
    gui.create_secondary_actors = False

    log = gui.log
    log.level = 'error'
    log.level = 'warning'
    #log.set_level('error')
    #log.set_level('warning')

    gui.load_nastran_geometry(bdf_filename)
    gui.load_nastran_results(op2_filename)
    vtk_ugrid = gui.grid
    if vtk_ugrid is None:
        raise RuntimeError('vtk_ugrid is None')

    save_nastran_results(gui, vtk_ugrid)
    #root = vtkMultiBlockDataSet()
    #coords_branch = vtkMultiBlockDataSet()

    #root.SetBlock(0, vtk_ugrid)
    #root.SetBlock(1, coords_branch)

    #for coord_id, axes_actor in test.axes.items():
        #coords_branch.SetBlock(0, axes)

    writer = vtkXMLUnstructuredGridWriter()
    writer.SetFileName(vtk_filename)
    writer.SetInputData(vtk_ugrid)
    out = writer.Write()
    #print('done', out)


def add_vtk_array(location: str,
                  point_data: vtkPointData,
                  cell_data: vtkCellData,
                  vtk_array) -> None:
    if location == 'node':
        point_data.AddArray(vtk_array)
    else:
        assert location == 'centroid'
        cell_data.AddArray(vtk_array)
    return

def main() -> None:  # pragma: no cover
    PKG_PATH = pyNastran.__path__[0]
    MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

    bdf_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.bdf')
    op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
    vtk_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.vtu')
    nastran_to_vtk(bdf_filename, op2_filename, vtk_filename)

    op2_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.op2')
    vtk_filename = os.path.join(MODEL_PATH, 'elements', 'modes_elements.vtu')
    nastran_to_vtk(op2_filename, op2_filename, vtk_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
    print("done")
