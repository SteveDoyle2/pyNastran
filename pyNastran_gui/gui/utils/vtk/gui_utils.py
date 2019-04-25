"""
these are functions not exposed in the GUI that are still useful
"""
def add_actors_to_gui(gui, actors, render=True):
    """adds multiple vtk actors"""
    if not len(actors):
        return
    renderer = gui.rend
    for actor in actors:
        renderer.AddActor(actor)
    if render:
        renderer.Render()

def remove_actors_from_gui(gui, actors, render=True):
    """removes multiple vtk actors"""
    if not len(actors):
        return
    renderer = gui.rend
    for actor in actors:
        renderer.RemoveActor(actor)
    if render:
        renderer.Render()
