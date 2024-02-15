import vtk

meta_image_filename = 'labels.mhd'

# Prepare to read the file

reader_volume = vtk.vtkMetaImageReader()
reader_volume.SetFileName(meta_image_filename)
reader_volume.Update()


# Extract the region of interest
voi = vtk.vtkExtractVOI()
if vtk.VTK_MAJOR_VERSION <= 5:
    voi.SetInput(reader_volume.GetOutput())
else:
    voi.SetInputConnection(reader_volume.GetOutputPort())

#voi.SetVOI(0,517, 0,228, 0,392)
voi.SetSampleRate(1, 1, 1)
#voi.SetSampleRate(3, 3, 3)
voi.Update() # necessary for GetScalarRange()
srange = voi.GetOutput().GetScalarRange() # needs Update() before!
print("Range", srange)


##Prepare surface generation
#contour = vtk.vtkContourFilter()
#contour = vtk.vtkMarchingCubes()
contour = vtk.vtkDiscreteMarchingCubes() #for label images
if vtk.VTK_MAJOR_VERSION <= 5:
    contour.SetInput(voi.GetOutput())
else:
    contour.SetInputConnection(voi.GetOutputPort())
contour.ComputeNormalsOn()


##run through all labels
for index in range(1, int(srange[1]) + 1):
    print("Doing label", index)

    contour.SetValue(0, index)
    contour.Update() #needed for GetNumberOfPolys() !!!

    smoother = vtk.vtkWindowedSincPolyDataFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        smoother.SetInput(contour.GetOutput())
    else:
        smoother.SetInputConnection(contour.GetOutputPort())
    smoother.SetNumberOfIterations(5)
    #smoother.BoundarySmoothingOff()
    #smoother.FeatureEdgeSmoothingOff()
    #smoother.SetFeatureAngle(120.0)
    #smoother.SetPassBand(.001)
    smoother.NonManifoldSmoothingOn()
    smoother.NormalizeCoordinatesOn()
    smoother.Update()


    ##calc cell normal
    triangle_cell_normals = vtk.vtkPolyDataNormals()
    if vtk.VTK_MAJOR_VERSION <= 5:
        triangle_cell_normals.SetInput(smoother.GetOutput())
    else:
        triangle_cell_normals.SetInputConnection(smoother.GetOutputPort())
    triangle_cell_normals.ComputeCellNormalsOn()
    triangle_cell_normals.ComputePointNormalsOff()
    triangle_cell_normals.ConsistencyOn()
    triangle_cell_normals.AutoOrientNormalsOn()
    triangle_cell_normals.Update() #creates vtkPolyData


    ##calc cell area
    triangle_cell_mesh_quality = vtk.vtkMeshQuality()
    if vtk.VTK_MAJOR_VERSION <= 5:
        triangle_cell_mesh_quality.SetInput(triangle_cell_normals.GetOutput())
    else:
        triangle_cell_mesh_quality.SetInputConnection(triangle_cell_normals.GetOutputPort())
    triangle_cell_mesh_quality.SetTriangleQualityMeasureToArea()
    triangle_cell_mesh_quality.SaveCellQualityOn() #default
    triangle_cell_mesh_quality.Update() #creates vtkDataSet

    point_normal_array = triangle_cell_normals.GetOutput().GetCellData().GetNormals()
    quality_array = triangle_cell_mesh_quality.GetOutput().GetCellData().GetArray("Quality")

    if point_normal_array.GetNumberOfTuples() != quality_array.GetNumberOfTuples():
        print("Error! Sizes of normal array and area array dont equal!")
        exit(1)

    f = open('label_stat' + "_%.4d" % index + ".dat", 'w')
    f.write("#cell_index\tarea\tn_x\tn_y\tn_z")

    for i in range(0, point_normal_array.GetNumberOfTuples()):

        point_normal = point_normal_array.GetTuple3(i) #this is for 3D data in python
        area = quality_array.GetValue(i)
        f.write('%s %s %s %s %s\n' % (i, area, point_normal[0], point_normal[1], point_normal[2]))

    f.close()
