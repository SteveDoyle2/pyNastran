import vtk
from pyNastran.gui.vtk_rendering_core import vtkRenderer, vtkRenderWindow, vtkActor, vtkPolyDataMapper

class vtkTimerCallback:
    def __init__(self):
        self.timer_count = 0

    def execute(self,obj,event):
        #print(self.timer_count)
        self.actor.SetPosition(self.timer_count, self.timer_count,0)
        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1


def main():
    #Create a sphere
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetCenter(0.0, 0.0, 0.0)
    sphere_source.SetRadius(5)

    #Create a mapper and actor
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(sphere_source.GetOutputPort())
    actor = vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()

    # Setup a renderer, render window, and interactor
    renderer = vtkRenderer()
    render_window = vtkRenderWindow()
    #render_window.SetWindowName('Test')

    render_window.AddRenderer(renderer)
    render_window_interactor = vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    #Add the actor to the scene
    renderer.AddActor(actor)
    renderer.SetBackground(1,1,1) # Background color white

    #Render and interact
    render_window.Render()

    # Initialize must be called prior to creating timer events.
    render_window_interactor.Initialize()

    # Sign up to receive TimerEvent
    cb = vtkTimerCallback()
    cb.actor = actor
    render_window_interactor.AddObserver('TimerEvent', cb.execute)
    timerId = render_window_interactor.CreateRepeatingTimer(100)

    #start the interaction and timer
    render_window_interactor.Start()


if __name__ == '__main__':  # pragma: no cover
    main()
