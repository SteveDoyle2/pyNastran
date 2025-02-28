"""
a transcription of:
 - https://blog.kitware.com/developing-hdf5-readers-using-vtkpythonalgorithm/
 - https://blog.kitware.com/vtkpythonalgorithm-is-great/
 - http://berkgeveci.github.io/2015/01/19/block-streaming/

 https://calcul.math.cnrs.fr/attachments/spip/IMG/pdf/vtk_visualization_pipeline.pdf
 https://cvw.cac.cornell.edu/ParaViewAdv/progsource

"""
from __future__ import annotations
import os
from itertools import chain
from collections import defaultdict

import numpy as np
import h5py
import vtk
import vtkmodules

from pyNastran.gui.vtk_interface import vtkUnstructuredGrid
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util import keys
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from pyNastran.dev.h5.read_h5 import pyNastranH5
from pyNastran.utils import object_methods, object_stats, object_attributes
from pyNastran.gui.utils.vtk.base_utils import numpy_to_vtk, numpy_to_vtkIdTypeArray
from pyNastran.dev.h5.fill_unstructured_grid import (
    fill_gui_vtk_unstructured_grid,
    numpy_to_vtk_points)

vtkInformation = vtkmodules.vtkCommonCore.vtkInformation
vtkPoints = vtkmodules.vtkCommonCore.vtkPoints
vtkQuad = vtkmodules.vtkCommonDataModel.vtkQuad
vtkInformationVector = vtkmodules.vtkCommonCore.vtkInformationVector


def write_file(data, xfreq, dirname: str):
    wdata = dsa.WrapDataObject(data)
    array = wdata.PointData['RTData']
    # Note that we flip the dimensions here because
    # VTK's order is Fortran whereas h5py writes in
    # C order. We don't want to do deep copies so we write
    # with dimensions flipped and pretend the array is
    # C order.
    array = array.reshape(wdata.GetDimensions()[::-1])
    hdf5_filename = os.path.join(dirname, 'data%d.h5' % xfreq)
    f = h5py.File(hdf5_filename, 'w')
    f.create_dataset("RTData", data=array)


def setup(dirname):
    rt = vtk.vtkRTAnalyticSource()
    for xfreq in range(60, 80):
        rt.SetXFreq(xfreq)
        rt.Update()
        write_file(rt.GetOutput(), xfreq)


metaDataKey = keys.MakeKey(
    keys.DataObjectMetaDataKey, "nastran_poly_data", "my module")
class HDF5Source(VTKPythonAlgorithmBase):
    """
    https://blog.kitware.com/a-vtk-pipeline-primer-part-1/
    https://blog.kitware.com/a-vtk-pipeline-primer-part-2/
    """
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=0,
            nOutputPorts=1, outputType='vtkPolyData')
        self.__FileName = ''

    def RequestInformation(self, request, inInfo, outInfo):
        print("MySource RequestInformation:")
        outInfo.GetInformationObject(0).Set(metaDataKey, vtk.vtkPolyData())
        print(outInfo.GetInformationObject(0))
        return 1
    def RequestData(self, request, inInfo, outInfo):
        info = outInfo.GetInformationObject(0)
        output = vtk.vtkPolyData.GetData(info)
        print(info)
        return 1

    def SetFileName(self, fname: str):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname
            assert isinstance(fname, str)

    def GetFileName(self) -> str:
        return self.__FileName

class HDF5Source(VTKPythonAlgorithmBase):
    """
    FillInputPortInformation(self, port, info) : Same as described above. Except self == vtkself.
    FillOutputPortInformation(self, port, info) : Same as described above. Except self == vtkself.
    RequestDataObject(self, request, inInfo, outInfo) : This is where you can create output data objects if the output DATA_TYPE_NAME() is not a concrete type.
    RequestInformation(self, request, inInfo, outInfo) : Provide meta-data downstream. More on this on later blogs.
    RequestUpdateExtent(self, request, inInfo, outInfo) : Modify requests coming from downstream. More on this on later blogs.
    RequestData(self, request, inInfo, outInfo) : Produce data. As described before.
    """
    def __init__(self):
        """
        1 input port + 1 output port == your common VTK filter.
        0 input port + 1 output port == your common VTK source.
        1 input port + 0 output port == your common VTK sink (writer, mapper etc.).
        """
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=0,
            nOutputPorts=1,
            outputType='vtkUnstructuredGrid')
        self.__FileName = ''

    #def FillOutputPortInformation(self, port, info):
        #info.Set(vtk.vtkDataObject.DATA_TYPE_NAME(), 'vtkUnstructuredGrid')
        #return 1

    def RequestData(self, request, inInfo, outInfo) -> int:
        filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\utils\hdf5\data60.h5'
        vtk_grid = vtkUnstructuredGrid()

        #f = h5py.File(filename, 'r')
        #data = f['RTData'][:]
        output = dsa.WrapDataObject(vtk_grid.GetData(outInfo))
        # Note that we flip the dimensions here because
        # VTK's order is Fortran whereas h5py writes in
        # C order.
        #output.SetDimensions(data.shape[::-1])

        #output.PointData.append(data.ravel(), 'RTData')
        #output.PointData.SetActiveScalars('RTData')
        output.PointData.append(nids, 'NodeIDs')
        output.PointData.SetActiveScalars('NodeIDs')

        vtk_info = outInfo.GetInformationObject(0)
        #vtk.info.Set(
            #vtk.vtkAlgorithm.CAN_PRODUCE_SUB_EXTENT(), 0)
        pipe = vtk.vtkCompositeDataPipeline
        #vtk_info.Set(
            #pipe.REQUEST_DATA_OBJECT(),
            #vtk_grid)
        vtk_info.Set(
            pipe.REQUEST_DATA(),
            vtk_grid)
        self._vtk_grid = vtk_grid
        #super().ProcessRequest(request, inputVector,outputVector)
        return 1

    def SetFileName(self, fname: str):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname
            assert isinstance(fname, str)

    def GetFileName(self) -> str:
        return self.__FileName

    def RequestInformation(self, request: vtkInformation,
                           inInfo,
                           outInfo: vtkInformationVector) -> int:
        """
        As I discussed previously, RequestInformation provides meta-data downstream.
        This meta-data is most of the time lightweight. In this example, we used
        f[‘RTData’].shape to read extent meta-data from the HDF5 file. This does not
        read any heavyweight data. Later on, we will see other examples of meta-data
        that is provided during RequestInformation.
        """

        self.model = pyNastranH5()
        self.model.read_h5_nastran(self.__FileName)

        filename = r'C:\NASA\m4\formats\git\pyNastran\pyNastran\utils\hdf5\data60.h5'
        f = h5py.File(filename, 'r')
        # Note that we flip the shape because VTK is Fortran order
        # whereas h5py reads in C order. When writing we pretend that the
        # data was C order so we have to flip the extents/dimensions.
        dims = f['RTData'].shape[::-1]  # (21, 21, 21)
        vtk_info = outInfo.GetInformationObject(0)

        #dims2 = (0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1)
        #vtk_info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(), dims2, 6)

            #vtk_info.Set(
                #vtk.vtkCompositeDataPipeline.REQUEST_DATA_OBJECT(),
                #vtk_grid)

        return 1

def fill_gui_vtk_unstructured_grid_results(model: pyNastranH5,
                                           vtk_ugrid: vtkUnstructuredGrid,
                                           eids: np.ndarray) -> None:
    for i, res in model.results.items():
        if res.location == 'node':
            names, resultsi = res.get_results()
            for name, result in zip(names, resultsi):
                print(name)

                # remove
                #point_data.SetActiveVectors(name)
                #point_data.SetVectors(result)
        elif res.location == 'element':
            names, resultsi = res.get_results()
            for name, result in zip(names, resultsi):
                print(name)

def set_caero_grid(alt_grids, ncaeros_points: int, model: BDF):
    """
    Sets the CAERO panel geometry.

    Parameters
    ----------
    ncaeros_points : int
        number of points used by the 'caero' actor
    model : BDF()
        the bdf model

    """
    log = model.log
    #gui = self.gui
    #vtk_points = vtk.vtkPoints()
    #vtk_points.SetNumberOfPoints(ncaeros_points)

    max_cpoints = []
    min_cpoints = []

    zfighting_offset = 0.0001
    caero_grid = alt_grids['caero']
    j = 0
    points = []
    for unused_eid, element in sorted(model.caeros.items()):
        if element.type in ('CAERO1', 'CAERO3', 'CAERO4', 'CAERO5', 'CAERO7'):
            # wing panel
            p1 = element.p1.copy()
            p4 = element.p4.copy()
            p1[0] += zfighting_offset
            p4[0] += zfighting_offset

            p2 = p1 + np.array([element.x12, 0., 0.])
            p3 = p4 + np.array([element.x43, 0., 0.])
            #[p1, p2, p3, p4]
            cpoints = [p1, p2, p3, p4]
            #cpoints = element.get_points()
            #cpoints[0][2] += zfighting_offset
            #cpoints[1][2] += zfighting_offset
            max_cpoints.append(np.array(cpoints).max(axis=0))
            min_cpoints.append(np.array(cpoints).min(axis=0))

            elem = vtkQuad()
            elem.GetPointIds().SetId(0, j)
            elem.GetPointIds().SetId(1, j + 1)
            elem.GetPointIds().SetId(2, j + 2)
            elem.GetPointIds().SetId(3, j + 3)

            points.extend(cpoints)
            #points.InsertPoint(j, *cpoints[0])
            #points.InsertPoint(j + 1, *cpoints[1])
            #points.InsertPoint(j + 2, *cpoints[2])
            #points.InsertPoint(j + 3, *cpoints[3])
            caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            j += 4
        elif element.type in ('CAERO2', 'BODY7'):
            # slender body

            #if 0:  # pragma: no cover
            # 1D version
            #cpoints = element.get_points()
            #cpoints[:, 2] += zfighting_offset
            #max_cpoints.append(np.array(cpoints).max(axis=0))
            #min_cpoints.append(np.array(cpoints).min(axis=0))

            #elem = vtk.vtkLine()
            #elem.GetPointIds().SetId(0, j)
            #elem.GetPointIds().SetId(1, j + 1)
            #points.InsertPoint(j, *cpoints[0])
            #points.InsertPoint(j + 1, *cpoints[1])
            #j += 2
            #caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            #else:
            # 3D version
            xyz, elems = element.get_points_elements_3d()
            assert xyz is not None, element
            xyz[:, 2] += zfighting_offset
            for elemi in elems:
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j + 1)
                elem.GetPointIds().SetId(2, j + 2)
                elem.GetPointIds().SetId(3, j + 3)
                n1, n2, n3, n4 = elemi

                points.append(xyz[n1])
                points.append(xyz[n2])
                points.append(xyz[n3])
                points.append(xyz[n4])
                #points.InsertPoint(j, *xyz[n1])
                #points.InsertPoint(j + 1, *xyz[n2])
                #points.InsertPoint(j + 2, *xyz[n3])
                #points.InsertPoint(j + 3, *xyz[n4])

                #cpoints = element.get_points()
                #cpoints[0][2] += zfighting_offset
                #cpoints[1][2] += zfighting_offset
                #max_cpoints.append(np.array(cpoints).max(axis=0))
                #min_cpoints.append(np.array(cpoints).min(axis=0))

                caero_grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j += 4
        else:
            log.info("skipping %s" % element.type)

    if ncaeros_points and len(max_cpoints):
        log.info('CAERO.max = %s' % np.vstack(max_cpoints).max(axis=0))
        log.info('CAERO.min = %s' % np.vstack(min_cpoints).min(axis=0))


    nodes = np.array(points, dtype='float32')
    vtk_points = numpy_to_vtk_points(nodes, points=None, dtype='<f', deep=1)
    caero_grid.SetPoints(vtk_points)

def set_spc_mpc_suport_grid(model: BDF,
                            alt_grids: dict[str, vtkUnstructuredGrid],
                            nid_to_pid_map: dict[int, int],
                            nid_map: dict[int, int],
                            idtype: str):
    """
    for each subcase, make secondary actors including:
     - spc_id=spc_id
     - mpc_id=mpc_id              (includes rigid elements)
     - mpc_dependent_id=mpc_id    (includes rigid elements)
     - mpc_independent_id=mpc_id  (includes rigid elements)
     - suport_id=suport1_id       (includes SUPORT/SUPORT1)

    TODO: consider changing the varying ids to huh???
    """
    spc_names = []
    mpc_names = []
    suport_names = []

    #print('getting rigid')
    rigid_lines = model._get_rigid()

    spc_ids_used = set()
    mpc_ids_used = set()
    suport1_ids_used = set()
    spc_to_subcase = defaultdict(list)
    mpc_to_subcase = defaultdict(list)
    #suport1_to_subcase = defaultdict(list)
    if 0:  # pragma: no cover
        for subcase_id, subcase in sorted(model.subcases.items()):
            if 'SPC' in subcase:
                spc_id = subcase.get_int_parameter('SPC')
                if spc_id is not None:
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0

                    ## TODO: this line seems too loose...
                    ## TODO: why aren't SPCDs included?
                    if nspcs + nspc1s + nspcds:
                        spc_to_subcase[spc_id].append(subcase_id)

            if 'MPC' in subcase:
                mpc_id = subcase.get_parameter('MPC')[0]
                if mpc_id is not None:

                    ## TODO: this line seems too loose
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs:
                        mpc_to_subcase[mpc_id].append(subcase_id)

        for spc_id in chain(model.spcs, model.spcadds):
            spc_name = f'SPC={spc_id:d}'
            if spc_id in mpc_to_subcase:
                subcases = spc_to_subcase[spc_id]
                spc_name += ': Subcases='
                spc_name += ', '.join(str(subcase_id) for subcase_id in subcases)
            spc_names += self._fill_spc(spc_id, spc_name, model, nid_to_pid_map)

        if nastran_settings.is_rbe:
            for mpc_id in chain(model.mpcs, model.mpcadds):
                depname = 'MPC=%i_dependent' % mpc_id
                indname = 'MPC=%i_independent' % mpc_id
                linename = 'MPC=%i_lines' % mpc_id
                if mpc_id in mpc_to_subcase:
                    subcases = mpc_to_subcase[mpc_id]
                    mpc_name = ': Subcases='
                    mpc_name += ', '.join(str(subcase_id) for subcase_id in subcases)
                    depname += mpc_name
                    indname += mpc_name
                    linename += mpc_name

            lines = get_mpc_node_ids(model, mpc_id, stop_on_failure=False)
            lines2 = list(lines)
            mpc_names += self._fill_dependent_independent(
                mpc_id, model, lines2,
                depname, indname, linename, idtype)

    if 0:  # pragma: no cover
        for subcase_id, subcase in sorted(model.subcases.items()):
            if 'SPC' in subcase:
                spc_id = subcase.get_int_parameter('SPC')
                if spc_id is not None and spc_id not in spc_ids_used:
                    spc_ids_used.add(spc_id)
                    nspcs = model.card_count['SPC'] if 'SPC' in model.card_count else 0
                    nspc1s = model.card_count['SPC1'] if 'SPC1' in model.card_count else 0
                    nspcds = model.card_count['SPCD'] if 'SPCD' in model.card_count else 0

                    ## TODO: this line seems too loose...
                    ## TODO: why aren't SPCDs included?
                    if nspcs + nspc1s + nspcds:
                        spc_name = 'spc_id=%i' % spc_id
                        spc_names += self._fill_spc(spc_id, spc_name, model, nid_to_pid_map)

            # rigid body elements and MPCs
            if 'MPC' in subcase:
                mpc_id = subcase.get_int_parameter('MPC')
                if mpc_id is not None and mpc_id not in mpc_ids_used:
                    mpc_ids_used.add(mpc_id)

                    ## TODO: this line seems too loose
                    nmpcs = model.card_count['MPC'] if 'MPC' in model.card_count else 0
                    if nmpcs and is_rbe:
                        lines = get_mpc_node_ids(model, mpc_id, stop_on_failure=False)
                        lines2 = list(lines)
                        depname = 'mpc_id=%i_dependent' % mpc_id
                        indname = 'mpc_id=%i_independent' % mpc_id
                        linename = 'mpc_id=%i_lines' % mpc_id
                        mpc_names += self._fill_dependent_independent(
                            mpc_id, model, lines2,
                            depname, indname, linename, idtype)

            # SUPORTs are node/dofs that deconstrained to allow rigid body motion
            # SUPORT1s are subcase-specific SUPORT cards
            if 'SUPORT1' in subcase.params:  ## TODO: should this be SUPORT?
                suport_id = subcase.get_int_parameter('SUPORT1')

                # TODO: is this line correct???
                if 'SUPORT' in model.card_count or 'SUPORT1' in model.card_count:

                    # TODO: this "if block" seems unnecessary
                    if suport_id is not None and suport_id not in suport1_ids_used:
                        # SUPORT1 / SUPORT
                        suport1_ids_used.add(suport_id)
                        suport_name = self._fill_suport(suport_id, subcase_id, model)
                        suport_names.append(suport_name)


    # create a SUPORT actor if there are no SUPORT1s
    # otherwise, we already included it in suport_id=suport_id
    if len(suport_names) == 0 and model.suport:
        # handle SUPORT without SUPORT1
        ids = []
        for suport in model.suport:
            idsi = suport.node_ids
            ids += idsi
        grid_name = 'SUPORT'
        suport_ugrid = vtkUnstructuredGrid()
        #i = np.searchsorted(model.nids, ids)
        #alt_grids[grid_name]
        #numpy_to_vtk_points(idsi)
        alt_grids[grid_name] = suport_ugrid
        #self.gui.create_alternate_vtk_grid(
            #grid_name, color=RED_FLOAT, opacity=1.0, point_size=4,
            #representation='point', is_visible=True)

    if len(rigid_lines):
        # handle RBEs without MPCs
        mpc_id = 0
        depname = 'rigid_dependent'
        indname = 'rigid_independent'
        linename = 'rigid_lines'
        mpc_names += _fill_dependent_independent(
            mpc_id, model, rigid_lines,
            depname, indname, linename, idtype,
            alt_grids, nid_map)

    geometry_names = spc_names + mpc_names + suport_names
    return geometry_names

def _fill_dependent_independent(unused_mpc_id: int, model: BDF,
                                lines,
                                depname: str, indname: str, linename: str,
                                idtype: str,
                                alt_grids: dict[str, vtkUnstructuredGrid],
                                nid_map: dict[int, int]) -> list[str]:
    """creates the mpc actors"""
    if not lines:
        return []

    #self.gui.create_alternate_vtk_grid(
        #depname, color=GREEN_FLOAT, line_width=5, opacity=1.,
        #point_size=5, representation='point', is_visible=False)
    #self.gui.create_alternate_vtk_grid(
        #indname, color=LIGHT_GREEN_FLOAT, line_width=5, opacity=1.,
        #point_size=5, representation='point', is_visible=False)
    #self.gui.create_alternate_vtk_grid(
        #linename, color=LIGHT_GREEN_FLOAT, line_width=5, opacity=1.,
        #point_size=5, representation='wire', is_visible=False)
    line_grid = vtkUnstructuredGrid()
    alt_grids[linename] = line_grid

    lines2 = []
    for line in lines:
        if line not in lines2:
            lines2.append(line)
    lines = np.array(lines2, dtype=idtype)
    dependent = (lines[:, 0])
    independent = np.unique(lines[:, 1])
    #self.dependents_nodes.update(dependent)
    unused_node_ids = np.unique(lines.ravel())

    #msg = ', which is required by %r' % depname
    #self._add_nastran_nodes_to_grid(depname, dependent, model, msg)

    #msg = ', which is required by %r' % indname
    #self._add_nastran_nodes_to_grid(indname, independent, model, msg)

    msg = ', which is required by %r' % linename
    _add_nastran_lines_to_grid(line_grid, nid_map,
                               linename, lines, model)

    mpc_names = [depname, indname, linename]
    return mpc_names

def _add_nastran_lines_to_grid(alt_grid: vtkUnstructuredGrid,
                               nid_map: dict[int, int],
                               name: str, lines, model: BDF, nid_to_pid_map=None):
    """used to create MPC lines"""
    nlines = lines.shape[0]
    #nids = np.unique(lines)
    #nnodes = len(nids)
    nnodes = nlines * 2
    if nnodes == 0:
        return
    #self.gui.follower_nodes[name] = lines.ravel()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nnodes)

    j = 0
    etype = 3 # vtkLine
    #nid_map = self.gui.nid_map
    #alt_grid = self.gui.alt_grids[name]
    for nid1, nid2 in lines:
        try:
            unused_i1 = nid_map[nid1]
        except KeyError:
            model.log.warning('nid=%s does not exist' % nid1)
            continue
        try:
            unused_i2 = nid_map[nid2]
        except KeyError:
            model.log.warning('nid=%s does not exist' % nid2)
            continue

        if nid1 not in model.nodes or nid2 not in model.nodes:
            continue
        node = model.nodes[nid1]
        point = node.get_position()
        points.InsertPoint(j, *point)

        node = model.nodes[nid2]
        point = node.get_position()
        points.InsertPoint(j + 1, *point)

        elem = vtk.vtkLine()
        point_ids = elem.GetPointIds()
        point_ids.SetId(0, j)
        point_ids.SetId(1, j + 1)
        alt_grid.InsertNextCell(etype, point_ids)
        j += 2
    alt_grid.SetPoints(points)

def fill_vtk_unstructured_grid_aero(model: BDF,
                                    nodes, node_ids,
                                    alt_grids: dict[str, vtkUnstructuredGrid]):
    if model.caeros:
        ncaeros_points = 0
        ugrid_aero = vtkUnstructuredGrid()
        alt_grids['caero'] = ugrid_aero
        set_caero_grid(alt_grids, ncaeros_points, model)

    if model.masses:
        set_mass_grid(alt_grids, model, nodes, node_ids)

def set_mass_grid(alt_grids: dict[str, vtkUnstructuredGrid],
                  model: BDF, nodes, node_ids) -> None:
    #nmass = len(model.masses)
    #eid_array = np.full(nmass, -1, dtype='int32')
    #mass_array = np.full(nmass, np.nan, dtype='float32')
    #xyz_array = np.full((nmass, 3), np.nan, dtype='float32')
    i = 0
    eids = []
    mass = []
    xyz = []
    mass_grid = vtkUnstructuredGrid()
    for eid, element in model.masses.items():
        if element.type == 'CONM2':
            nid = element.nid
            inid = np.searchsorted(node_ids, nid)
            xyz_nid = nodes[inid, :]
            centroid = element.offset(xyz_nid)

            #eid_array[i] = eid
            #mass_array[i] = element.mass
            #xyz_array[i, :] = centroid
            eids.append(eid)
            mass.append(element.mass)
            xyz.append(centroid)
            vtk_elem = vtk.vtkVertex()
            vtk_elem.GetPointIds().SetId(0, inid)
            mass_grid.InsertNextCell(vtk_elem.GetCellType(), vtk_elem.GetPointIds())
        else:
            model.log.debug(f'skipping:\n{str(elem)}')
            continue
        i += 1
    eid_array = np.array(eids, dtype='int32')[:i]
    mass_array = np.array(mass, dtype='float32')[:i]
    xyz_array = np.array(xyz, dtype='float32')[:i, :]

    vtk_points = numpy_to_vtk_points(xyz_array, points=None, dtype='<f', deep=1)
    mass_grid.SetPoints(vtk_points)

    # add vectors
    vtk_eids = numpy_to_vtk(eid_array, deep=0, array_type=None)
    vtk_mass = numpy_to_vtk(mass_array, deep=0, array_type=None)
    vtk_eids.SetName('Mass ID')
    vtk_mass.SetName('Mass')
    #point_data = mass_grid.GetPointData()
    cell_data = mass_grid.GetCellData()
    cell_data.AddArray(vtk_mass)
    alt_grids['mass'] = mass_grid

def fill_vtk_unstructured_grid_constraints(model: BDF,
                                           alt_grids: dict[str, vtkUnstructuredGrid],
                                           nid_map: dict[int, int]):
    #if model.caeros:
        #ncaeros_points = 0
        #ugrid_aero = vtkUnstructuredGrid()
        #alt_grids['caero'] = ugrid_aero
    nid_to_pid_map = None
    idtype = 'int32'
    set_spc_mpc_suport_grid(model,
                            alt_grids,
                            nid_to_pid_map,
                            nid_map,
                            idtype)

def get_gui_nastran_ugrid(hdf5_filename: str,
                          ugrid_main: vtkUnstructuredGrid,
                          add_property_info: bool=True,
                          add_material_info: bool=True,
                          subcases=None,  # default=None -> all
                          modes=None, # default=None -> all
                          results=None, # default=None -> all,
                          ) -> tuple[BDF, vtkUnstructuredGrid]:
    add_aero = False
    add_constraints = False
    add_results = False

    model = pyNastranH5(add_aero, add_constraints, add_results, subcases)
    model.read_h5_nastran(hdf5_filename)
    geom_model = model.geom_model
    geom_model.log.info(geom_model.card_count)

    #ugrid_main = vtkUnstructuredGrid()
    alt_grids = {
        'main' : ugrid_main,
    }
    node_ids, eids, form, cases = fill_gui_vtk_unstructured_grid(
        geom_model, ugrid_main,
        add_property=add_property_info,
        add_material=add_material_info,
    )
    #fill_vtk_unstructured_grid_aero(geom_model, nodes, node_ids, alt_grids)
    #fill_vtk_unstructured_grid_constraints(geom_model, alt_grids, nid_map)
    fill_gui_vtk_unstructured_grid_results(model, ugrid_main, eids)


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
    return model, ugrid_main, root, alt_grids, node_ids, eids, form, cases

def add_actor_to_renderer(actor: vtk.vtkActor):
    ren = vtk.vtkRenderer()
    ren.AddActor(actor)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(ren)
    render_window.SetSize(600, 600)
    render_window.Render()

    # Setup render window interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    #style = vtk.vtkInteractorStyleImage()

    # Render and start interaction
    render_window_interactor.SetRenderWindow(render_window)
    render_window_interactor.Initialize()

    render_window_interactor.Start()
    return ren, render_window, render_window_interactor

    #for xfreq in range(60, 80):
        #hdf5_filename = os.path.join(dirname, 'data%d.h5' % xfreq)
        #alg.SetFileName(hdf5_filename)
        #renWin.Render()
    #asdf

    #import vtk, h5py
    #from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
    #from vtk.numpy_interface import dataset_adapter as dsa
