import os
import wx
import vtk
from vtk.wx.wxVTKRenderWindow import wxVTKRenderWindow
from vtk import (vtkTriangle,vtkQuad,vtkTetra,vtkWedge,vtkHexahedron,
                 vtkQuadraticTriangle,vtkQuadraticQuad,vtkQuadraticTetra,vtkQuadraticWedge,vtkQuadraticHexahedron)
from numpy import zeros,ones

import pyNastran
version = pyNastran.__version__

from pyNastran.bdf.bdf import *
from pyNastran.op2.op2 import OP2
from mouseStyle import MouseStyle


def getScreenCorner(x,y):
    #print "wx.GetDisplaySize() = ",wx.GetDisplaySize()
    (xScreen,yScreen) = wx.GetDisplaySize()
    xCorner = (xScreen-x)//2
    yCorner = (yScreen-y)//2
    return(xCorner,yCorner)

class pyWidget(wxVTKRenderWindow):
    def __init__(self,*args,**kwargs):
        wxVTKRenderWindow.__init__(self, *args, **kwargs)
        self.parent = args[0]
        self.dirname = ""
        self.OnChar = self.onChar2

    def ResetCamera(self):
        self.Reset()

    def GetCamera(self):
        return self._CurrentCamera

    def onChar2(self,event):
        #print "onChar2 = ",event.GetKeyCode()
        camera = self.GetCamera()
        code = event.GetKeyCode()
        if   code == ord('m'): # zooming in
            camera.Zoom(1.1)
        elif code == ord('M'): # zooming out
            camera.Zoom(0.9)

        elif code == ord('o'): # counter-clockwise
            camera.Roll(5.)
        elif code == ord('O'): # clockwise
            camera.Roll(-5.)

        elif code == ord('x'): # set x-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(    0.,0., 1.)
            camera.SetPosition(  1.,0., 0.)
            self.ResetCamera()
        elif code == ord('X'): # set x-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(    0.,0.,-1.)
            camera.SetPosition( -1.,0., 0.)
            self.ResetCamera()


        elif code == ord('y'): # set y-axis
            camera.SetFocalPoint(0.,0.,0.)
            camera.SetViewUp(    0.,0.,1.)
            camera.SetPosition(  0.,1.,0.)
            self.ResetCamera()
        elif code == ord('Y'): # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(    0., 0.,-1.)
            camera.SetPosition(  0.,-1., 0.)
            self.ResetCamera()

        elif code == ord('z'): # set z-axis
            camera.SetFocalPoint(0.,0.,0.)
            camera.SetViewUp(    0.,1.,0.)
            camera.SetPosition(  0.,0.,1.)
            self.ResetCamera()
        elif code == ord('Z'): # set z-axis
            camera.SetFocalPoint(0.,0., 0.)
            camera.SetViewUp(   0., -1.,0.)
            camera.SetPosition( 0., 0.,-1.)
            self.ResetCamera()

        elif code == ord('i'):
            self.TakePicture(event)

        elif code == ord('e'): # edges dont work right yet
            self.parent.DisplayEdges(event)

        elif code == ord('L'):
            self.parent.cycleResults()

        self.Update()
        self.Render()
        ###

    def TakePicture(self,event):
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.getRenderer())
        renderLarge.SetMagnification(4)

        wildcard = "PNG (*.png)|*.png|" \
         "JPEG (*.jpeg; *.jpeg; *.jpg; *.jfif)|*.jpg;*.jpeg;*.jpe;*.jfif|" \
         "TIFF (*.tif; *.tiff)|*.tif;*.tiff|" \
         "BMP (*.bmp)|*.bmp|" \
         "PostScript (*.ps)|*.ps|" \
         "All files (*.*)|*.*"
        
        dlg = wx.FileDialog(None, "Choose a file", self.dirname, "", wildcard, wx.SAVE | wx.OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            fname        = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            fname = os.path.join(self.dirname,fname)

            print "fname = ",fname

            # We write out the image which causes the rendering to occur. If you
            # watch your screen you might see the pieces being rendered right
            # after one another.
            lfname = fname.lower()
            if lfname.endswith('.png'):
                writer = vtk.vtkPNGWriter()
            elif lfname.endswith('.jpeg'):
                writer = vtk.vtkJPEGWriter()
            elif lfname.endswith('.tiff'):
                writer = vtk.vtkTIFFWriter()
            elif lfname.endswith('.ps'):
                writer = vtk.vtkPostScriptWriter()
            else:
                writer = vtk.vtkPNGWriter()

            writer.SetInputConnection(renderLarge.GetOutputPort())
            writer.SetFileName(fname)
            writer.Write()
        dlg.Destroy()

    def getRenderer(self):
        return self.GetCurrentRenderer()


class Pan(wx.Panel):
    def __init__(self, *args, **kwargs):
        wx.Panel.__init__(self, *args, **kwargs)
        isEdges = False
        self.isEdges = isEdges # surface wireframe
        self.widget = pyWidget(self, -1)

        window = self.widget.GetRenderWindow()
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(window)

    def DisplayEdges(self,event):
        self.isEdges = not(self.isEdges)
        if 0:
            self.getEdges()

        if 1:
            prop = self.edgeActor.GetProperty()
            print "dir(prop) = ",dir(prop)
            print "visible = ",prop.GetEdgeVisibility()
            if self.isEdges:
                prop.EdgeVisibilityOn()
                print "edges are now on\n"
            else:
                prop.EdgeVisibilityOff()
                print "edges are now off\n"
            prop.Modified()
            #prop.Update()
            #self.edgeActor.Modified()
            #self.edgeActor.Update()
        if 0:
            print "dir(widget) = ",dir(self.widget)

            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                print "iactor = ",i
                actor = actors.GetNextItem()
                
                try:
                    if self.isEdges:
                        actor.getProperty().EdgeVisibilityOn()
                    else:
                        actor.getProperty().EdgeVisibilityOff()
                    print "set the edges..."
                except:
                    pass
                #actor.GetProperty().SetLineWidth(0.0)
        self.widget.Render()
        self.widget.Update()


    def SetToWireframe(self,event):
        if self.bdfFileName is not None:
            self.widget.Wireframe()
            self.widget.Update()

    def SetToSurface(self,event):
        if self.bdfFileName is not None:
            self.widget.Surface()
            self.widget.Update()

    def SetToFlatShading(self,event): # Flat
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToFlat()
            self.widget.Render()
        ###

    def SetToGouraudShading(self,event): # Gouraud
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToGouraud()
            self.widget.Render()
        ###

    def SetToPhongShading(self,event): # Phong
        if self.bdfFileName is not None:
            actors = self.widget._CurrentRenderer.GetActors()
            numActors = actors.GetNumberOfItems()
            actors.InitTraversal()
            for i in range(0,numActors):
                actor = actors.GetNextItem()
                actor.GetProperty().SetInterpolationToPhong()
            self.widget.Render()
        ###

    def WireframeTemp(self):
        """Sets the current actor representation as wireframe.
        """
        actors = self.widget._CurrentRenderer.GetActors()
        numActors = actors.GetNumberOfItems()
        actors.InitTraversal()
        for i in range(0,numActors):
            actor = actors.GetNextItem()
            actor.GetProperty().SetRepresentationToWireframe()
        self.Render()

    def getColors(self):
        pass
        #Merger = vtk.vtkMergeFilter()
        #Filter = vtk.vtkGeometryFilter()
        #Filter.SetInput(self.aQuadGrid)
        
        #Filter = vtk.vtkRotationFilter()
        #Filter.SetCenter(0.,0.,0.)
        #Filter.SetNumberOfCopies(1)
        #Filter.SetInput(aQuadGrid)
        #Filter.Update()
        
    def getEdges(self):
        edges = vtk.vtkExtractEdges()
        edges.SetInput(self.grid)
        self.edgeMapper = vtk.vtkPolyDataMapper()
        self.edgeMapper.SetInput(edges.GetOutput())
        #self.edgeMapper.EdgeVisibilityOff()

        self.edgeActor = vtk.vtkActor()
        self.edgeActor.SetMapper(self.edgeMapper)
        self.edgeActor.GetProperty().SetColor(0,0,0)
        prop = self.edgeActor.GetProperty()
        #prop.SetLineWidth(0.0)
        if self.isEdges:
            prop.EdgeVisibilityOn()
        else:
            prop.EdgeVisibilityOff()

        self.rend.AddActor(self.edgeActor)

        print "visible = ",prop.GetEdgeVisibility()
        #self.edgeActor.Update()


    def addGeometry(self):
        print "addGeometry"
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid)
        
        #lut = vtk.vtkLookupTable()
        #aQuadMapper.SetLookupTable(lut)

        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 0, 0) # red
        #geometryActor.GetProperty().SetBackfaceProperty(1, 0, 0) # red
        #geometryActor.GetProperty().BackfaceCullingOn()  # hidges elements that have normals not facing camera
        #geometryActor.GetProperty().SetLineWidth(0.5)
        self.rend.AddActor(geometryActor)

    def addAltGeometry(self):
        print "addAltGeometry"
        aQuadMapper = vtk.vtkDataSetMapper()
        aQuadMapper.SetInput(self.grid2)
        #aQuadMapper.SetInput(Filter.GetOutput())
        geometryActor = vtk.vtkActor()
        geometryActor.SetMapper(aQuadMapper)
        #geometryActor.AddPosition(2, 0, 2)
        #geometryActor.GetProperty().SetDiffuseColor(0, 0, 1) # blue
        geometryActor.GetProperty().SetDiffuseColor(1, 1, 0) # green
        self.rend.AddActor(geometryActor)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def Update(self):
        #print '\n'.join(dir(self.widget))
        window = self.widget.GetRenderWindow()
        self.rend.ResetCamera()
        window.Render()
        self.widget.Update()
        #self.rend.Update()
        #self.renWin.Render()
        
    def main(self):
        print "main builder"
        window = self.widget.GetRenderWindow()
        self.rend = vtk.vtkRenderer()
        window.AddRenderer(self.rend)

        self.addGeometry()
        self.addAltGeometry()

        # Create the usual rendering stuff.
        if self.isEdges:
            self.getEdges()
        self.rend.GetActiveCamera().ParallelProjectionOn()
        self.rend.SetBackground(.1,.2,.4)
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        self.rend.ResetCamera()

        #iren.Start()

    def getWindow(self):
        return self.widget.GetRenderWindow()

    def getWindowName(self,winName=''):
        return "pyNastran v%s - %s" %(version,self.bdfFileName)
        #return "pyNastran v%s - %s" %(version,'solid.bdf')

    def setWindowName(self,winName=''):
        window = self.getWindow()
        window.SetWindowName("pyNastran v%s - %s" %(version,self.bdfFileName))
    
    def isStress(self,op2,ID):
        if ID in op2.solidStress or ID in op2.plateStress or ID in op2.compositePlateStress or ID in op2.barStress  or ID in op2.beamStress or ID in op2.rodStress:
            return True
        return False
        
    def loadResults(self,op2FileName):
        op2 = OP2(op2FileName)
        op2.readOP2()
        #print op2.printResults()
        
        #case = op2.displacements[1]
        #print "case = ",case
        #for nodeID,translation in sorted(case.translations.items()):
        #    print "nodeID=%s t=%s" %(nodeID,translation)
        #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]

        cases = {}
        subcaseIDs = op2.iSubcaseNameMap.keys()
        self.iSubcaseNameMap = op2.iSubcaseNameMap
        
        eidsSet = False
        for ID in subcaseIDs:
            if not eidsSet:
                eids = zeros(self.nElements,'d')
                for eid,eid2 in self.eidMap.items():
                    eids[eid2] = float(eid)
               
                eKey = (ID,'isElementOn',1,'centroid')
                cases[(ID,'Element_ID',1,'centroid')] = eids
                cases[eKey] = ones(self.nElements) # is the element supported
                eidsSet = True

            if ID in op2.displacements:
                case = op2.displacements[ID]
                key = (ID,'DisplacementX',3,'node')
                cases[key] = case.translations

            if self.isStress(op2,ID):
                cases = self.fillStressCase(cases,op2,ID,eKey)
            ###
            
        ###
        self.resultCases = cases
        self.caseKeys = sorted(cases.keys())
        print "caseKeys = ",self.caseKeys
        #print "type(caseKeys) = ",type(self.caseKeys)
        self.iCase = -1
        self.nCases = len(self.resultCases)-1 # number of keys in dictionary
        self.cycleResults() # start at nCase=0

    def fillStressCase(self,cases,op2,ID,eKey):
        oxx = zeros(self.nElements)
        oyy = zeros(self.nElements)
        ozz = zeros(self.nElements)

        o1  = zeros(self.nElements)
        o2  = zeros(self.nElements)
        o3  = zeros(self.nElements)
        ovm = zeros(self.nElements)

        if ID in op2.rodStress:
            case = op2.rodStress[ID]
            for eid in case.axial:
                print "bar eid=%s" %(eid)
                axial   = case.axial[eid]
                torsion = case.torsion[eid]

                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 0.
                oxx[eid2] = axial
                oyy[eid2] = torsion

                o1[eid2]  = max(axial,torsion)  # not really
                #oyy[eid2] = torsion
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = min(axial,torsion)

        if ID in op2.barStress:
            #self.s1    = {}
            #self.s2    = {}
            #self.s3    = {}
            #self.s4    = {}
            case = op2.barStress[ID]
            for eid in case.axial:
                print "bar eid=%s" %(eid)
                oxxi = case.axial[eid]
                o1i  = max(case.smax[eid])
                o3i  = min(case.smin[eid])
                
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 0.

                oxx[eid2] = oxxi
                #oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                #o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                #ovm[eid2] = ovmi

        if ID in op2.plateStress:
            #self.txy    = {}

            case = op2.plateStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                print "plate eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                #self.oyy[eid][nid][iLayer]

                oxxi = max(case.oxx[eid]['C'])
                oyyi = max(case.oyy[eid]['C'])
                ozzi = min(case.oxx[eid]['C'],min(case.oyy[eid]['C']))

                o1i = max(case.majorP[eid]['C'])
                o2i = max(case.minorP[eid]['C'])
                o3i = min(case.majorP[eid]['C'],min(case.minorP[eid]['C']))
                ovmi = max(case.ovmShear[eid]['C'])


                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 0.
                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                #ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = 0.  #(o1i+o3i)/2.
                o3[eid2] = o3i
                ovm[eid2] = ovmi

        if ID in op2.solidStress:
            case = op2.solidStress[ID]
            if case.isVonMises():
                vmWord = 'vonMises'
            else:
                vmWord = 'maxShear'
            for eid in case.ovmShear:
                print "solid eid=%s" %(eid)
                oxxi = case.oxx[eid]['C']
                oyyi = case.oyy[eid]['C']
                ozzi = case.ozz[eid]['C']
                o1i = case.o1[eid]['C']
                o2i = case.o2[eid]['C']
                o3i = case.o3[eid]['C']
                ovmi = case.ovmShear[eid]['C']
                #if ID==1:
                #    print "ovm[%s] = %s" %(eid,ovmi)
                eid2 = self.eidMap[eid]
                cases[eKey][eid2] = 0.
                oxx[eid2] = oxxi
                oyy[eid2] = oyyi
                ozz[eid2] = ozzi

                o1[eid2] = o1i
                o2[eid2] = o2i
                o3[eid2] = o3i
                ovm[eid2] = ovmi
            ###

        cases[(ID,'StressXX',1,'centroid')] = oxx
        cases[(ID,'StressYY',1,'centroid')] = oyy
        cases[(ID,'StressZZ',1,'centroid')] = ozz

        cases[(ID,'Stress1',1,'centroid')] = o1
        cases[(ID,'Stress2',1,'centroid')] = o2
        cases[(ID,'Stress3',1,'centroid')] = o3
        cases[(ID,vmWord   ,1,'centroid')] = ovm
        return cases

    def incrementCycle(self):
        if self.iCase is not self.nCases:
            self.iCase +=1
        else:
            self.iCase = 0
        key = self.caseKeys[self.iCase]
        if key[2] ==3: # vector size=3 -> vector
            self.incrementCycle()
        
        #print "next key = ",key
        return

    def cycleResults(self):
        self.incrementCycle()
        self.gridResult = vtk.vtkFloatArray()

        #allocationSize = vectorSize*location (where location='node'-> self.nNodes)
        #self.gridResult.Allocate(self.nNodes,1000)

        #allocationSize = vectorSize*location (where location='centroid'-> self.nElements)
        self.gridResult.Allocate(self.nElements,1000)

        key = self.caseKeys[self.iCase]
        case = self.resultCases[key]

        (subcaseID,resultType,vectorSize,location) = key

        #self.iSubcaseNameMap[self.iSubcase] = [Subtitle,Label]
        caseName = self.iSubcaseNameMap[subcaseID]
        subtitle,label = caseName

        print "subcaseID=%s resultType=%s subtitle=%s label=%s" %(subcaseID,resultType,subtitle,label)

        for value in case:
            maxValue = value
            minValue = value
            break

        for value in case:
            maxValue = max(value,maxValue)
            minValue = min(value,minValue)
        
        # flips sign to make colors go from blue -> red
        normValue = maxValue-minValue
        print "case = ",case
        if normValue==0.: # avoids division by 0.
            normValue = 1.
        for value in case:
            self.gridResult.InsertNextValue((value-minValue)/normValue)
        print "max=%g min=%g norm=%g\n" %(maxValue,minValue,normValue)

        self.grid.GetCellData().SetScalars(self.gridResult)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        self.grid.Modified()

    def loadGeometry(self,bdfFileName,dirname):

        if bdfFileName is None:
            self.grid       = vtk.vtkUnstructuredGrid()
            self.gridResult = vtk.vtkFloatArray()
            #self.vectorResult = vtk.vtkFloatArray()
            self.grid2      = vtk.vtkUnstructuredGrid()
            return
        else:
            self.grid.Reset()
            self.gridResult.Reset()
            
            self.grid2.Reset()
            try:
                del self.resultCases
                del self.caseKeys
                del self.iCase
                del self.nCases
                del self.iSubcaseNameMap
            except:
                print "cant cleanup..."
            ###

        model = BDF()
        model.readBDF(bdfFileName,includeDir=dirname)

        nNodes    = model.nNodes()
        nElements = model.nElements()
        nCAeros   = model.nCAeros()
        self.nNodes = nNodes
        self.nElements = nElements

        #print "nNodes = ",nNodes
        print "nElements = ",nElements

        #self.aQuadGrid.Allocate(nElements+nNodes, 1000)

        if 'CONM2' in model.cardCount:
            nCONM2 = model.cardCount['CONM2']
        else:
            nCONM2 = 0
        self.grid.Allocate(nElements, 1000)
        self.grid2.Allocate(nCAeros +nCONM2, 1000)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(nNodes)
        self.gridResult.Allocate(nNodes,1000)
        #vectorReselt.SetNumberOfComponents(3)
        self.nidMap = {}
        i=0
        #elem.SetNumberOfPoints(nNodes)
        fraction = 1./nNodes
        for nid,node in sorted(model.nodes.items()):
            #print "i = ",i
            point = node.Position()
            #print "point = ",point
            #sys.stdout.flush()
            #aVoxel = vtk.vtkPixel()
            #print "made voxel"; sys.stdout.flush()
            points.InsertPoint(i, *point)
            self.gridResult.InsertNextValue(i*fraction)
            #print str(element)

            #elem = vtk.vtkVertex()
            #elem.GetPointIds().SetId(0, i)
            #self.aQuadGrid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            #vectorResult.InsertTuple3(0, 0.0, 0.0, 1.0)

            self.nidMap[nid] = i
            i+=1

        j = 0
        points2 = vtk.vtkPoints()
        points2.SetNumberOfPoints(nCAeros*4+nCONM2)
        for eid,element in sorted(model.caeros.items()):
            if isinstance(element,CAERO1):
                cpoints = element.Points()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, j)
                elem.GetPointIds().SetId(1, j+1)
                elem.GetPointIds().SetId(2, j+2)
                elem.GetPointIds().SetId(3, j+3)
                points2.InsertPoint(j,   *cpoints[0])
                points2.InsertPoint(j+1, *cpoints[1])
                points2.InsertPoint(j+2, *cpoints[2])
                points2.InsertPoint(j+3, *cpoints[3])
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=4
        self.mapElements(points,points2,self.nidMap,model,j)

    def mapElements(self,points,points2,nidMap,model,j):
        self.eidMap = {}
        i = 0
        for eid,element in sorted(model.elements.items()):
            self.eidMap[eid] = i
            if isinstance(element,CTRIA3):
                #print "ctria3"
                elem = vtkTriangle()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTRIA6):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTriangle()
                    elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                else:
                    elem = vtkTriangle()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD4):
                nodeIDs = element.nodeIDs()
                elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CQUAD8):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticQuad()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                else:
                    elem = vtkQuad()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA4):
                elem = vtkTetra()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CTETRA10):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticTetra()
                    elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                    elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                    elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8, nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9, nidMap[nodeIDs[9]])
                else:
                    elem = vtkTetra()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CPENTA6):
                elem = vtkWedge()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())

            elif isinstance(element,CPENTA15):
                nodeIDs = element.nodeIDs()
                if None not in nodeIDs:
                    elem = vtkQuadraticWedge()
                    elem.GetPointIds().SetId(6,  nidMap[nodeIDs[6]])
                    elem.GetPointIds().SetId(7,  nidMap[nodeIDs[7]])
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                else:
                    elem = vtkWedge()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA8):
                nodeIDs = element.nodeIDs()
                elem = vtkHexahedron()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,CHEXA20):
                nodeIDs = element.nodeIDs()
                #print "nodeIDs = ",nodeIDs
                if None not in nodeIDs:
                    elem = vtkQuadraticHexahedron()
                    elem.GetPointIds().SetId(8,  nidMap[nodeIDs[8]])
                    elem.GetPointIds().SetId(9,  nidMap[nodeIDs[9]])
                    elem.GetPointIds().SetId(10, nidMap[nodeIDs[10]])
                    elem.GetPointIds().SetId(11, nidMap[nodeIDs[11]])
                    elem.GetPointIds().SetId(12, nidMap[nodeIDs[12]])
                    elem.GetPointIds().SetId(13, nidMap[nodeIDs[13]])
                    elem.GetPointIds().SetId(14, nidMap[nodeIDs[14]])
                    elem.GetPointIds().SetId(15, nidMap[nodeIDs[15]])
                    elem.GetPointIds().SetId(16, nidMap[nodeIDs[16]])
                    elem.GetPointIds().SetId(17, nidMap[nodeIDs[17]])
                    elem.GetPointIds().SetId(18, nidMap[nodeIDs[18]])
                    elem.GetPointIds().SetId(19, nidMap[nodeIDs[19]])
                else:
                    elem = vtkHexahedron()

                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                elem.GetPointIds().SetId(2, nidMap[nodeIDs[2]])
                elem.GetPointIds().SetId(3, nidMap[nodeIDs[3]])
                elem.GetPointIds().SetId(4, nidMap[nodeIDs[4]])
                elem.GetPointIds().SetId(5, nidMap[nodeIDs[5]])
                elem.GetPointIds().SetId(6, nidMap[nodeIDs[6]])
                elem.GetPointIds().SetId(7, nidMap[nodeIDs[7]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            elif isinstance(element,LineElement) or isinstance(element,SpringElement):
                elem = vtk.vtkLine()
                nodeIDs = element.nodeIDs()
                elem.GetPointIds().SetId(0, nidMap[nodeIDs[0]])
                elem.GetPointIds().SetId(1, nidMap[nodeIDs[1]])
                self.grid.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
            ###
            elif isinstance(element,CONM2): # not perfectly located
                nid  = element.Nid()
                c    = element.Centroid()
                elem = vtk.vtkVertex()
                #elem = vtk.vtkSphere()
                #elem.SetRadius(1.0)
                #print str(element)

                points2.InsertPoint(j,     *c)
                elem.GetPointIds().SetId(0, j)
                #elem.SetCenter(points.GetPoint(nidMap[nid]))
                self.grid2.InsertNextCell(elem.GetCellType(), elem.GetPointIds())
                j+=1
            else:
                print "skipping %s" %(element.type)
            i+=1
        ###
        self.grid.SetPoints(points)
        self.grid2.SetPoints(points2)
        #self.grid.GetPointData().SetScalars(self.gridResult)
        #self.grid.GetCellData().SetScalars(self.gridResult)
        self.grid.Modified()
        self.grid2.Modified()
        self.grid.Update()
        self.grid2.Update()
        print "updated grid"

    def OnKeyPress(self,obj,event):
        rwi = obj
        key = rwi.GetKeySym()
        print "*Pressed %s" %(key)

    def buildVTK(self,bdfFileName=None,dirname=None):
        self.bdfFileName = bdfFileName
        self.loadGeometry(self.bdfFileName,dirname)
        self.main()

        #self.renWin = vtk.vtkRenderWindow()
        #renWin = window
        
        cam = self.rend.GetActiveCamera()
        mouseArgs = {'pipeline':self,'camera':cam}

        xSize = 800
        ySize = 600
        (x,y) = getScreenCorner(xSize,ySize)
        #window.SetSize(xSize,ySize)
        #window.SetPosition(x,y)

        #iren = wxVTKRenderWindowInteractor(self,-1)
        iren = vtk.vtkRenderWindowInteractor()
        mouseStyle = MouseStyle(mouseArgs,iren)
        iren.SetInteractorStyle(mouseStyle)
        window = self.widget.GetRenderWindow()
        iren.SetRenderWindow(window)

        iren.AddObserver("KeyPressEvent", self.OnKeyPress)

        #print "type(ren) = ",type(self.rend)
        #self.rend.GetActiveCamera().Zoom(2.0)

        #self.setWindowName()
        # Render the scene and start interaction.
        iren.Initialize()
        #iren.Start()
        #window.Render()

