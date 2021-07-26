"""tests the NastranIO class"""
import os
from copy import deepcopy
from typing import Set
import unittest

import numpy as np
#try:
    #import matplotlib
    #matplotlib.use('Agg')
    #IS_MATPLOTLIB = True
#except ModuleNotFoundError:  # pyparsing is missing
    #IS_MATPLOTLIB = False
#except ImportError:
    #pass
import vtk
from cpylog import SimpleLogger

import pyNastran
from pyNastran.bdf.bdf import BDF
#from pyNastran.bdf.cards.test.test_aero import get_zona_model
#from pyNastran.bdf.errors import DuplicateIDsError
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


def nastran_to_vtk(op2_filename: str, vtk_filename: str):
    test = NastranGUI()
    test.create_secondary_actors = False

    log = test.log
    log.level = 'error'

    test.load_nastran_geometry(op2_filename)
    test.load_nastran_results(op2_filename)
    vtk_ugrid = test.grid
    cases = test.result_cases.items()


    point_data = vtk_ugrid.GetPointData()
    cell_data = vtk_ugrid.GetCellData()

    used_titles = set()
    for key, case_data in cases:
        case, index_name = case_data
        if case.is_complex:
            log.warning(f'skipping case {str(case)} because it is complex')
            continue
        if isinstance(case, (GridPointForceResult)):
            log.warning(f'skipping case {str(case)}')
            continue

        #key
        #0

        #case_data[0]
        #GuiResult
            #title='NodeID'
            #data_type='<i4'
            #uname='GuiResult'

        #case_data[1]
        #(0, 'NodeID')


        #self.data_map = data_map
        #self.subcase_id = subcase_id
        ##assert self.subcase_id > 0, self.subcase_id

        #self.title = title
        #self.header = header
        ##self.scale = scale
        #self.location = location
        #assert location in ['node', 'centroid'], location
        #self.subcase_id = subcase_id
        #self.uname = uname

        #self.scalar = scalar
        ##self.data_type = self.dxyz.dtype.str # '<c8', '<f4'
        #self.data_type = self.scalar.dtype.str # '<c8', '<f4'
        #self.is_real = True if self.data_type in REAL_TYPES else False
        #self.is_complex = not self.is_real
        #self.nlabels = nlabels
        #self.labelsize = labelsize
        #self.ncolors = ncolors
        #self.colormap = colormap

        ##print('title=%r data_type=%r' % (self.title, self.data_type))
        #if self.data_type in INT_TYPES:
            #self.data_format = '%i'
        #elif data_format is None:
            #self.data_format = '%.2f'
        #else:
            #self.data_format = data_format

        #self.title_default = self.title
        #self.header_default = self.header
        #self.data_format_default = self.data_format

        #self.min_default = np.nanmin(self.scalar)
        #self.max_default = np.nanmax(self.scalar)
        #if self.data_type in INT_TYPES:
            ## turns out you can't have a NaN/inf with an integer array
            ## we need to recast it
            #if mask_value is not None:
                #inan_short = np.where(self.scalar == mask_value)[0]
                #if len(inan_short):
                    ## overly complicated way to allow us to use ~inan to invert the array
                    #inan = np.in1d(np.arange(len(self.scalar)), inan_short)
                    #inan_remaining = self.scalar[~inan]

                    #self.scalar = np.asarray(self.scalar, 'f')
                    #self.data_type = self.scalar.dtype.str
                    #self.data_format = '%.0f'
                    #self.scalar[inan] = np.nan
                    #try:
                        #self.min_default = inan_remaining.min()
                    #except ValueError:  # pragma: no cover
                        #print('inan_remaining =', inan_remaining)
                        #raise
                    #self.max_default = inan_remaining.max()
        #else:
            ## handling VTK NaN oddinty
            ## filtering the inf values and replacing them with NaN
            ## 1.#R = inf
            ## 1.#J = nan
            #ifinite = np.isfinite(self.scalar)
            #if not np.all(ifinite):
                #self.scalar[~ifinite] = np.nan
                #try:
                    #self.min_default = self.scalar[ifinite].min()
                #except ValueError:
                    #print(self.title)
                    #print(self.scalar)
                    #raise
                #self.max_default = self.scalar[ifinite].max()
        #self.min_value = self.min_default
        #self.max_value = self.max_default

        if isinstance(case, ForceTableResults):
            assert len(case.titles) == 1, case
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
                    dxyz = case.dxyz[itime, :, :]
                    vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
                    titlei =  f'{title}_subcase={case.subcase_id}'
                    _check_title(title, used_titles)
                    vtk_array.SetName(titlei)
                    _add_array(case.location, point_data, cell_data, vtk_array)
            else:
                log.warning(f'cannot add {str(case)}')
            continue
        elif isinstance(case, DisplacementResults):
            assert len(case.titles) == 1, case
            assert case.dxyz.ndim == 3, case.dxyz.shape
            for itime, title in enumerate(case.titles):
                dxyz = case.dxyz[itime, :, :]
                assert dxyz.ndim == 2, dxyz.shape
                vtk_array = numpy_to_vtk(dxyz, deep=0, array_type=None)
                titlei =  f'{title}_subcase={case.subcase_id}'
                _check_title(title, used_titles)
                vtk_array.SetName(titlei)
                _add_array(case.location, point_data, cell_data, vtk_array)
            continue
        elif isinstance(case, SimpleTableResults):
            log.warning(f'skipping SimpleTableResults {case}')
            continue
        elif isinstance(case, LayeredTableResults):
            log.warning(f'skipping LayeredTableResults {case}')
            continue
            assert case.location == 'centroid', case
            #(1, (0, 0, 'Static')) = index_name
            name = index_name[1]
            print(case, name)
            (itime, imethod, zero, unused_header) = name
            method = case.methods[imethod]
            #for method in case.methods:
            titlei =  f'{case.uname}_{method}_subcase={case.subcase_id}'                #assert titlei not in used_titles, titlei

            res = case.get_result(key, name)
            _check_title(titlei, used_titles)
            vtk_array = numpy_to_vtk(res, deep=0, array_type=None)
            vtk_array.SetName(titlei)
        else:
            title = case.title
            if title == 'Normals' and case.scalar is None:
                continue
            vtk_array = numpy_to_vtk(case.scalar, deep=0, array_type=None)
            if title in used_titles:
                log.warning(f'skipping {title} because it is already used')
                continue
            _check_title(title, used_titles)
            vtk_array.SetName(title)
        _add_array(case.location, point_data, cell_data, vtk_array)

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


def _check_title(title: str, used_titles: Set[str]) -> None:
    assert title not in used_titles, title
    used_titles.add(title)

def main():
    PKG_PATH = pyNastran.__path__[0]
    MODEL_PATH = os.path.join(PKG_PATH, '..', 'models')

    op2_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.op2')
    vtk_filename = os.path.join(MODEL_PATH, 'elements', 'static_elements.vtu')
    nastran_to_vtk(op2_filename, vtk_filename)

if __name__ == '__main__':  # pragma: no cover
    main()
    print("done")
