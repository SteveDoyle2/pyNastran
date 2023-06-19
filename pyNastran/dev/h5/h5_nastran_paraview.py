from __future__ import annotations
import os
from typing import TYPE_CHECKING

import numpy as np
import h5py
import vtk

from pyNastran.gui.vtk_renering_core import vtkDataSetMapper # , vtkPolyDataMapper
from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from pyNastran.dev.h5.fill_unstructured_grid import fill_paraview_vtk_unstructured_grid
from pyNastran.dev.h5.h5_nastran2 import add_actor_to_renderer, pyNastranH5
from pyNastran.dev.h5.vtk_request_subset import vtkRequestSubset
from pyNastran.dev.h5.vtk_h5_nastran import vtkH5NastranReader
if TYPE_CHECKING:  # pragma: no cover
    from pyNastran.bdf.bdf import BDF

def fill_paraview_vtk_unstructured_grid_results(model: pyNastranH5,
                                                vtk_ugrid: vtkUnstructuredGrid,
                                                eids: np.ndarray) -> None:
    point_data = vtk_ugrid.GetPointData()
    cell_data = vtk_ugrid.GetCellData()
    for i, res in model.results.items():
        if res.location == 'node':
            names, resultsi = res.get_results()
            for name, result in zip(names, resultsi):
                print(name)
                #name = 'cat'
                result.SetName(name)
                point_data.AddArray(result)

                # remove
                #point_data.SetActiveVectors(name)
                #point_data.SetVectors(result)
        elif res.location == 'element':
            res.eids = eids
            names, resultsi = res.get_results()
            for name, result in zip(names, resultsi):
                print(name)
                #name = 'cat'
                result.SetName(name)
                cell_data.AddArray(result)
                #cell_data.SetScalars(result)
        else:
            raise NotImplementedError(res)
        #break  # TODO: remove

def get_paraview_nastran_ugrid(hdf5_filename: str,
                               add_property_info=True,
                               add_material_info=True,
                               subcases=None,  # default=None -> all
                               modes=None, # default=None -> all
                               results=None, # default=None -> all,
                      ) -> tuple[BDF, vtkUnstructuredGrid]:
    #subcases = [1, 2]
    #modes = range(1, 10)
    #results = ['eigenvectors', 'displacement', 'stress', 'strain']
    add_aero = False
    add_constraints = False
    add_results = False

    model = pyNastranH5(add_aero, add_constraints, add_results, subcases)
    model.read_h5_nastran(hdf5_filename)
    geom_model = model.geom_model
    geom_model.log.info(geom_model.card_count)

    ugrid_main = vtkUnstructuredGrid()
    alt_grids = {
        'main' : ugrid_main,
    }
    eids = fill_paraview_vtk_unstructured_grid(
        geom_model, ugrid_main,
        add_property=add_property_info,
        add_material=add_material_info,
    )
    #fill_vtk_unstructured_grid_aero(geom_model, nodes, node_ids, alt_grids)
    #fill_vtk_unstructured_grid_constraints(geom_model, alt_grids, nid_map)
    fill_paraview_vtk_unstructured_grid_results(model, ugrid_main, eids)


    # PART 1 Make some Data.
    # Make a tree.
    root = vtk.vtkMultiBlockDataSet()

    #branch = vtk.vtkMultiBlockDataSet()
    #root.SetBlock(0, branch)

    # Make some leaves.
    #leaf1 = vtk.vtkSphereSource()
    #leaf1.SetCenter(0, 0, 0)
    #leaf1.Update()
    #branch.SetBlock(0, leaf1.GetOutput())
    #for key in ['mass']:
        #if key in alt_grids:
            #del alt_grids[key]

    iblock = 0
    basepath = os.path.basename(hdf5_filename)
    for name, ugrid in alt_grids.items():
        print(name)
        root.SetBlock(iblock, ugrid)
        meta_data = root.GetMetaData(iblock)
        meta_data.Set(vtk.vtkCompositeDataSet.NAME(), f'{basepath}: {name}')
        iblock += 1
    #root.SetBlock(1, grid_aero)

    #root.GetMetaData(0).Set(vtk.vtkCompositeDataSet.NAME(), basepath + ': main')
    #root.GetMetaData(1).Set(vtk.vtkCompositeDataSet.NAME(), basepath + ': CAERO')
    #print('root')
    return model, ugrid_main, root, alt_grids


def run_vtk(hdf5_filename: str, scale: float):
    #alg = HDF5Source()
    alg = vtkH5NastranReader()
    alg.SetFileName(hdf5_filename)

    #cf = vtk.vtkContourFilter()
    #cf.SetInputConnection(alg.GetOutputPort())
    #cf.SetValue(0, 200)
    model, ugrid, root, alt_grids = get_paraview_nastran_ugrid(hdf5_filename)
    ugrid = alt_grids['main']

    warp = vtk.vtkWarpVector()
    warp.SetScaleFactor(scale)
    warp.SetInputData(ugrid)
    warp.Update()

    #warp = ugrid
    grid_mapper = vtkDataSetMapper()
    if 0:
        grid_mapper.SetInputData(ugrid)
    else:
        grid_mapper.SetInputData(warp.GetOutput())

    #grid_mapper = vtkPolyDataMapper()
    #grid_mapper.SetInputConnection(ugrid.GetOutputPort())

    actor = vtk.vtkLODActor()
    actor.SetMapper(grid_mapper)

    add_actor_to_renderer(actor)
    print('go')
    #import time
    #while 1:
        #time.sleep(0.1)

def run_vtk2(hdf5_filename: str):
    #alg = HDF5Source()
    alg = vtkH5NastranReader()
    alg.SetFileName(hdf5_filename)

    rs = vtkRequestSubset()
    rs.SetInputConnection(alg.GetOutputPort())
    rs.SetUpdateExtent((5, 10, 5, 10, 0, 20))

    import vtk
    cf = vtk.vtkContourFilter()
    cf.SetInputConnection(rs.GetOutputPort())
    cf.SetValue(0, 200)

    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(cf.GetOutputPort())

    a = vtk.vtkActor()
    a.SetMapper(m)

    actor = a
    add_actor_to_renderer(actor)


def main():
    #hdf5_filename = r'C:\NASA\m4\formats\git\examples\car\GS61_PT_mode_TEST-noRB3.h5'
    #hdf5_filename = r'C:\NASA\m4\formats\git\examples\car\DM22_Transmission_mode_TEST.h5'
    #hdf5_filename = r'C:\NASA\m4\formats\git\examples\car\body-sol103_orig.h5'; scale = 5000.
    hdf5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\msc\mode_echo.h5'; scale = 1.0
    #hdf5_filename = r'C:\NASA\m4\formats\git\pyNastran\models\bwb\bwb_saero_saved.h5'; scale = 10.
    #setup(dirname)
    #run_vtk(hdf5_filename)
    import time
    t0 = time.time()
    run_vtk(hdf5_filename, scale)
    dt = time.time() - t0
    print('dt = ', dt)


if __name__ == '__main__':
    main()
