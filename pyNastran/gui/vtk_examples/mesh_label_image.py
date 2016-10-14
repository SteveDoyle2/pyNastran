import vtk

input='labels.mhd'

# Prepare to read the file

readerVolume = vtk.vtkMetaImageReader()
readerVolume.SetFileName(input)
readerVolume.Update()


# Extract the region of interest
voi = vtk.vtkExtractVOI()
if vtk.VTK_MAJOR_VERSION <= 5:
    voi.SetInput(readerVolume.GetOutput())
else:
    voi.SetInputConnection(readerVolume.GetOutputPort())

#voi.SetVOI(0,517, 0,228, 0,392)
voi.SetSampleRate(1,1,1)
#voi.SetSampleRate(3,3,3)
voi.Update()#necessary for GetScalarRange()
srange= voi.GetOutput().GetScalarRange()#needs Update() before!
print "Range", srange


##Prepare surface generation
#contour = vtk.vtkContourFilter()
#contour = vtk.vtkMarchingCubes()
contour = vtk.vtkDiscreteMarchingCubes() #for label images
if vtk.VTK_MAJOR_VERSION <= 5:
    contour.SetInput( voi.GetOutput() )
else:
    contour.SetInputConnection( voi.GetOutputPort() )
contour.ComputeNormalsOn()


##run through all labels

for index in range(1, int(srange[1]) + 1):

    print "Doing label", index

    contour.SetValue(0, index)
    contour.Update() #needed for GetNumberOfPolys() !!!


    smoother= vtk.vtkWindowedSincPolyDataFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        smoother.SetInput(contour.GetOutput());
    else:
        smoother.SetInputConnection(contour.GetOutputPort());
    smoother.SetNumberOfIterations(5);
    #smoother.BoundarySmoothingOff();
    #smoother.FeatureEdgeSmoothingOff();
    #smoother.SetFeatureAngle(120.0);
    #smoother.SetPassBand(.001);
    smoother.NonManifoldSmoothingOn();
    smoother.NormalizeCoordinatesOn();
    smoother.Update();


    ##calc cell normal
    triangleCellNormals= vtk.vtkPolyDataNormals()
    if vtk.VTK_MAJOR_VERSION <= 5:
        triangleCellNormals.SetInput(smoother.GetOutput())
    else:
        triangleCellNormals.SetInputConnection(smoother.GetOutputPort())
    triangleCellNormals.ComputeCellNormalsOn()
    triangleCellNormals.ComputePointNormalsOff()
    triangleCellNormals.ConsistencyOn()
    triangleCellNormals.AutoOrientNormalsOn()
    triangleCellNormals.Update() #creates vtkPolyData


    ##calc cell area
    triangleCellAN= vtk.vtkMeshQuality()
    if vtk.VTK_MAJOR_VERSION <= 5:
        triangleCellAN.SetInput(triangleCellNormals.GetOutput())
    else:
        triangleCellAN.SetInputConnection(triangleCellNormals.GetOutputPort())
    triangleCellAN.SetTriangleQualityMeasureToArea()
    triangleCellAN.SaveCellQualityOn() #default
    triangleCellAN.Update() #creates vtkDataSet

    PointNormalArray = triangleCellNormals.GetOutput().GetCellData().GetNormals()
    qualityArray = triangleCellAN.GetOutput().GetCellData().GetArray("Quality")

    if PointNormalArray.GetNumberOfTuples() != qualityArray.GetNumberOfTuples():
        print "Error! Sizes of normal array and area array dont equal!"
        exit(1)

    f= open('label_stat' + "_%.4d" % index + ".dat", 'w')
    print >> f, "#cell_index\tarea\tn_x\tn_y\tn_z"

    for i in range(0, PointNormalArray.GetNumberOfTuples()):

        pointNormal= PointNormalArray.GetTuple3(i) #this is for 3D data in python
        area= qualityArray.GetValue(i)
        print >> f, i, area, pointNormal[0], pointNormal[1], pointNormal[2]

    f.close()
