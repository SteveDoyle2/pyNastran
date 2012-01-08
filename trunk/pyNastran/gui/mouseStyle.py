import vtk

class MouseStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self,rwi,parent=None):
        self.AddObserver("MiddleButtonPressEvent",self.middleButtonPressEvent)
        self.AddObserver("MiddleButtonReleaseEvent",self.middleButtonReleaseEvent)
        self.interactor = rwi

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

    def middleButtonPressEvent(self,obj,event):
        print "Middle Button pressed"
        self.OnMiddleButtonDown()
        return

    def middleButtonReleaseEvent(self,obj,event):
        print "Middle Button released"
        self.OnMiddleButtonUp()
        return

    def OnKeyPress(self,asf):
       # Get the keypress
       rwi = self.Interactor
       key = rwi.GetKeySym()
  
       # Output the key that was pressed
       print "Pressed %s" %(key)
  
       # Handle an arrow key
       if(key == "Up"):
            print "Up"
 
       # Handle a "normal" key
       if(key == "a"):
           print "A"
 
       # Forward events
       self.OnKeyPress()
       rwi.Update()

