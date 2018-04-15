# kills the program when you hit Cntl+C from the command line
# doesn't save the current state as presumably there's been an error
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import random
import vtk
assert int(vtk.VTK_VERSION[0]) > 6, vtk.VTK_VERSION
#assert vtk.VTK_VERSION > (7, 1, 0), vtk.VTK_VERSION


def ActorCallback(caller, eventId_unuseD, clientData, callData_unused):

    text_actor = vtk.vtkBillboardTextActor3D(clientData)
    #actor = vtkActor *(caller)

    xyz = actor.GetPosition()
    label = '%.3f, %.3f, %.3f' % (xyz[0], xyz[1], xyz[2])
    text_actor.SetPosition(xyz)
    text_actor.SetInput(label)

def make_billboard(renderer, xyz):
    text_actor = vtk.vtkBillboardTextActor3D()
    #actor = vtkActor *(caller)
    #print(xyz)

    label = '%.3f, %.3f, %.3f' % xyz
    text_actor.SetPosition(xyz)
    text_actor.SetInput(label)

    #text_actor.SetPosition (actor.GetPosition())
    text_actor.GetTextProperty().SetFontSize(15)
    #text_actor.GetTextProperty().SetColor(1.0, 1.0, 0.4)
    text_actor.GetTextProperty().SetColor(0.0, 0.0, 0.0)
    text_actor.GetTextProperty().SetJustificationToCentered()

    renderer.AddActor(text_actor)

def main():
    # For testing
    #vtkMath::RandomSeed(8775070)

    # Create a sphere
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetCenter(0.0, 0.0, 0.0)
    sphere_source.SetRadius(1.0)

    # Create an actor
    mapper2 = vtk.vtkPolyDataMapper()
    mapper2.SetInputConnection(sphere_source.GetOutputPort())
    actor2 = vtk.vtkActor()
    actor2.SetMapper(mapper2)
    actor2.SetPosition(0, 0, 0)
    actor2.GetProperty().SetColor(1.0, .4, .4)

    # Create a renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(.6, .4, .2)
    renderer.AddActor(actor2)

    # Create a render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Create an interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    for i in range(10):
        # Create a mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere_source.GetOutputPort())

        # Create an actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.SetPosition(0, 0, 0)

        # Setup the text and add it to the renderer
        #text_actor = vtk.vtkBillboardTextActor3D()
        #text_actor.SetInput("")
        #text_actor.SetPosition (actor.GetPosition())
        #text_actor.GetTextProperty().SetFontSize (12)
        #text_actor.GetTextProperty().SetColor (1.0, 1.0, .4)
        #text_actor.GetTextProperty().SetJustificationToCentered()

        renderer.AddActor(actor)
        #renderer.AddActor(text_actor)

        actor_callback = vtk.vtkCallbackCommand()
        #actor_callback.SetCallback(ActorCallback)
        #actor_callback.SetClientData(textActor)
        #actor.AddObserver(vtkCommand::ModifiedEvent, actor_callback)
        xyz = (
            random.uniform(-10.0, 10.0),
            random.uniform(-10.0, 10.0),
            random.uniform(-10.0, 10.0)
        )
        actor.SetPosition(xyz)
        make_billboard(renderer, xyz)

    render_window.Render()
    render_window_interactor.Start()

if __name__ == '__main__':
    main()
