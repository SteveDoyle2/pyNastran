# Setup Note #
Download the entire package from subversion or just the [GUI](https://sourceforge.net/projects/pynastran/files/?source=navbar) executable.

If you download the source, make sure you follow the [Installation Guide](http://code.google.com/p/pynastran/wiki/InstallationGuide) and use **setup.py develop** and not **setup.py install**.

If you're using Python 3.x, the GUI will **not** work because [VTK](http://www.vtk.org/Wiki/VTK/Python_Wrapping_FAQ) does not support Python 3.

# Introduction #

The Graphical User Interface (GUI) looks like:

![https://pynastran.googlecode.com/svn/trunk/pyNastran/gui/qt.png](https://pynastran.googlecode.com/svn/trunk/pyNastran/gui/qt.png)

# Running the GUI #
On the command line:
```
>> pyNastranGUI

To view the options:
>> pyNastranGUI --help

Usage:
  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT]
                  [-s SHOT] [-m MAGNIFY] [-p SCRIPT]
                  [-q] [-e] [-n | -c]
  pyNastranGUI -h | --help
  pyNastranGUI -v | --version

Options:
  -h, --help                  show this help message and exit
  -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair,
                                           stl, tetgen, usm3d)
  -i INPUT, --input INPUT     path to input file
  -o OUTPUT, --output OUTPUT  path to output file
  -p SCRIPT, --pyscript SCIPRT  path to script file
  -s SHOT, --shots SHOT       path to screenshot (only 1 for now)
  -m MAGNIFY, --magnify       MAGNIFY how much should the resolution on a picture be magnified (default=1)
  -q, --quiet                 prints debug messages (default=True)
  -e, --edges                 shows element edges as black lines (default=False)
  -n, --nodalResults          plots nodal results (default)
  -c, --centroidalResults     plots centroidal results
  -v, --version               show program's version number and exit
```

The **solidBending.bdf** and **solidBending.op2** files have been included as examples that work in the GUI.  They are inside the "models" folder (at the same level as setup.py).

# Features #
  * Nastran, Panair, Cart3d, STL, Usm3d support
    * limited results support
  * Command line interface
  * Scripting capability
  * Take a Screenshot (menu/button/keyboard)
  * Snap to Axis (keyboard)
  * Change Background Color (menu)

# Supported Elements #
  * CQUAD4 / CQUAD8
  * CTRIA3 / CTRIA6
  * CTETRA4 / CTETRA10
  * CHEXA8 / CHEXA20
  * CPENTA6 / CPENTA15
  * CSHEAR
  * CQUADR / CTRIAR
  * CBAR / CBEAM / CROD / CONROD / CELASx (displayed as lines)
  * CAERO1 (shown in yellow)
  * CONM2 (shown in yellow as points)

# BDF Requirements #
  * Entire model can be cross-referenced
  * Same requirements as BDF (include an executive/case control deck, define all cross-referenced cards, etc.)

# Scripting #
GUI commands are logged to the window with their call signature.  Users may then use a custom Python script to take many pictures, show the sub-caero panels, etc.  A sample CAERO script that shows individual CAERO subpanels (instead of just the outline of the CAERO panel) is provided with the download.

## Graphical Issues ##
You'll have the best performance if you run the GUI on Windows with an new NVIDIA graphics card and on a desktop.

If you're having issues, you should update the driver for your graphics card, especially if you have a laptop or Radeon card. For a desktop machine, go to the web site of the manufacturer of the graphics card. For a laptop, you should normally go to the web site of the laptop manufacturer, though for NVIDIA you may now find a newer driver available from NVIDIA.

Issues include:

  1. the backfaces of elements not being colored
> 2. the GUI not working