## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
import os
import wx
import vtk
from vtk.wx.wxVTKRenderWindow import wxVTKRenderWindow

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

        #elif code == ord('e'): # edges dont work right yet
            #self.parent.DisplayEdges(event)

        elif code == ord('L'):
            self.parent.cycleResults()

        elif code == ord('h'):
            self.ShowHideScalarBar()

        self.Update()
        self.Render()
        ###

    def ShowHideScalarBar(self):
        if self.parent.nCases==0:
            return
        isOn = self.parent.scalarBar.GetVisibility()
        if isOn:
            self.parent.scalarBar.VisibilityOff()
            self.parent.TurnTextOff()
        else:
            self.parent.scalarBar.VisibilityOn()
            self.parent.TurnTextOn()
        self.parent.scalarBar.Modified()

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
