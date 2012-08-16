import vtk


class MouseStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, ren, parent=None):
        self.AddObserver("MiddleButtonPressEvent", self.middleButtonPressEvent)
        self.AddObserver("MiddleButtonReleaseEvent", self.middleButtonReleaseEvent)
        self.pipeline = ren['pipeline']
        self.cam = ren['camera']

    def Update(self):
        self.pipeline.widget.Modified()
        self.cam.Modified()
        self.cam.UpdateViewport(self.pipeline.rend)
        self.pipeline.renWin.Render()

    def getActiveCamera(self):
        ren = self.pipeline.rend
        #print "type(ren) = ",ren
        #print "dir = ",'\n'.join(dir(ren))
        camera = self.cam
        #camera = ren.getActiveCamera
        #print "type(camera) = ",camera
        return camera
        #return self.pipeline.rend.getActiveCamera()

    #def leftButtonPressEvent(self,obj,event):
    #    print "Left Button pressed"
    #    self.LeftButtonDown()
    #    return
    #
    #def rightButtonPressEvent(self,obj,event):
    #    print "Right Button pressed"
    #    self.OnRightButtonDown()
    #    return

    #def rightButtonReleaseEvent(self,obj,event):
    #    print "Right Button released"
    #    self.OnRightButtonUp()
    #    return

    def middleButtonPressEvent(self, obj, event):
        print "Middle Button pressed"
        self.OnMiddleButtonDown()
        return

    def middleButtonReleaseEvent(self, obj, event):
        print "Middle Button released"
        self.OnMiddleButtonUp()
        return

    def onChar(self, obj, event):
        rwi = obj
        key = rwi.GetKeySym()
        print "*Pressed %s" % (key)

        #renderer = self.ren
        camera = self.getActiveCamera()
        #print "type(camera) = ",type(camera)
        if key == 'm':  # zooming in
            camera.Zoom(1.1)
            self.Update()
        elif key == 'M':  # zooming out
            camera.Zoom(0.9)
            self.Update()
        elif key == 'd':  # draw edges
            if self.pipeline.isEdges == False:
                return
            prop = self.pipeline.edgeActor.GetProperty()
            #print '\n'.join(dir(prop))
            if self.pipeline.isEdges:
                #prop.SetLineWidth(0.0)
                prop.EdgeVisibilityOff()
            else:
                #prop.SetLineWidth(1.0)
                prop.EdgeVisibilityOn()
            #prop.Update()
            #prop.Modified()
            #prop.SetVisibility(False)
            self.pipeline.isEdges = not(self.pipeline.isEdges)
            self.pipeline.edgeMapper.Modified()
            self.pipeline.edgeMapper.Update()
            self.pipeline.edgeActor.Modified()

        # Roll
        elif key == 'o':  # counter-clockwise
            camera.Roll(5.)
            self.Update()
        elif key == 'O':  # clockwise
            camera.Roll(-5.)
            self.Update()

        # Yaw
        #elif key=='a': # counter-clockwise
            #camera.Yaw(5.)
            #self.Update()
        #elif key=='A': # clockwise
            #camera.Yaw(-5.)
            #self.Update()

        # Elevation
        #elif key=='v': # counter-clockwise
            #camera.Elevation(5.)
            #self.Update()
        #elif key=='V': # clockwise
            #camera.Elevation(-5.)
            #self.Update()

        # Pitch
        #elif key=='c': # counter-clockwise
            #camera.Pitch(5.)
            #self.Update()
        #elif key=='C': # clockwise
            #camera.Pitch(-5.)
            #self.Update()

        elif key == 'x':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(1., 0., 0.)
            self.pipeline.rend.ResetCamera()
            self.Update()
        elif key == 'X':  # set x-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(-1., 0., 0.)
            self.pipeline.rend.ResetCamera()
            self.Update()

        elif key == 'y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., 1.)
            camera.SetPosition(0., 1., 0.)
            self.pipeline.rend.ResetCamera()
            self.Update()
        elif key == 'Y':  # set y-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 0., -1.)
            camera.SetPosition(0., -1., 0.)
            self.pipeline.rend.ResetCamera()
            self.Update()

        elif key == 'z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., 1., 0.)
            camera.SetPosition(0., 0., 1.)
            self.pipeline.rend.ResetCamera()
            self.Update()
        elif key == 'Z':  # set z-axis
            camera.SetFocalPoint(0., 0., 0.)
            camera.SetViewUp(0., -1., 0.)
            camera.SetPosition(0., 0., -1.)
            self.pipeline.rend.ResetCamera()
            self.Update()

        #elif key=='i': # picture taking doesnt work
            #self.pipeline.takePicture()

        # Panning doesnt work
        #elif key=='Up':
            #p = camera.GetPosition()
            #print "p = ",p
            #f = camera.GetFocalPoint()
            #print "f = ",f
            #camera.SetFocalPoint()
            #camera.SetPosition(p[0],p[1],p[2]+1.)
            #camera.SetFocalPoint(f[0],f[1],f[2]+1.)
            #camera.Pan()
            #self.Update()
        #elif key=='Down':
            #p = camera.GetPosition()
            #print "p = ",p
            #f = camera.GetFocalPoint()
            #print "f = ",f
            #camera.SetFocalPoint()
            #camera.SetPosition(p[0],p[1],p[2]-1.)
            #camera.SetFocalPoint(f[0],f[1],f[2]-1.)
            #camera.Pan()
            #self.pipeline.rend.ResetCamera()
            #self.Update()
