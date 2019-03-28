"""
these are functions not exposed in the GUI that are still useful
"""
def add_actors(gui, actors, render=True):
    """adds multiple vtk actors"""
    renderer = gui.rend
    for actor in actors:
        renderer.AddActor(actor)
    if render:
        renderer.Render()

def remove_actors(gui, actors, render=True):
    """removes multiple vtk actors"""
    renderer = gui.rend
    for actor in actors:
        renderer.RemoveActor(actor)
    if render:
        renderer.Render()
