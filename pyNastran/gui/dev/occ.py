from OCC.Display.SimpleGui import init_display
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox

display, start_display, add_meu, add_function_to_menu = init_display()
mybox = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()
display.DisplayShape(mybox, update=True)
start_display()
