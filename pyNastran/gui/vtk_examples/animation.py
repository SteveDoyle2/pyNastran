from __future__ import print_function
import vtk

class vtkTimerCallback():
    def __init__(self):
        self.timer_count = 0

    def execute(self,obj,event):
        print(self.timer_count)
        self.actor.SetPosition(self.timer_count, self.timer_count,0);
        iren = obj
        iren.GetRenderWindow().Render()
        self.timer_count += 1


def main():
    #Create a sphere
    sphereSource = vtk.vtkSphereSource()
    sphereSource.SetCenter(0.0, 0.0, 0.0)
    sphereSource.SetRadius(5)

    #Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphereSource.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()

    # Setup a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    #renderWindow.SetWindowName("Test")

    renderWindow.AddRenderer(renderer);
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    #Add the actor to the scene
    renderer.AddActor(actor)
    renderer.SetBackground(1,1,1) # Background color white

    #Render and interact
    renderWindow.Render()

    # Initialize must be called prior to creating timer events.
    renderWindowInteractor.Initialize()

    # Sign up to receive TimerEvent
    cb = vtkTimerCallback()
    cb.actor = actor
    renderWindowInteractor.AddObserver('TimerEvent', cb.execute)
    timerId = renderWindowInteractor.CreateRepeatingTimer(100);

    #start the interaction and timer
    renderWindowInteractor.Start()


if __name__ == '__main__':
    main()