Graphical User Interface (GUI)
==============================

Setup Note
----------
Download the entire package from Github or just the `GUI
<https://sourceforge.net/projects/pynastran/files/?source=navbar/>`_ executable.

If you download the source, make sure you follow the `Installation Guide
<https://github.com/SteveDoyle2/pyNastran/wiki/Installation>`_ and use
**setup.py develop** and not **setup.py install**.

Python 2.7 with ``vtk==5.10.1`` will probably give you the best looking GUI.
VTK 6 and 7 probably still have some issues.
The GUI in Python 3 won't save your settings.

Introduction
------------

The Graphical User Interface (GUI) looks like:

.. image:: ../../pyNastran/gui/qt.png

The GUI also has a sidebar and transient support.

Running the GUI
---------------
On the command line:

.. code-block:: console

  >>> pyNastranGUI

To view the options:

.. code-block:: console

  >>> pyNastranGUI --help

    pyNastranGUI [-f FORMAT] INPUT
                 [-s SHOT] [-m MAGNIFY]
                 [-g GSCRIPT] [-p PSCRIPT]
                 [-u POINTS_FNAME...]
                 [-q]
    pyNastranGUI [-f FORMAT] INPUT OUTPUT
                 [-s SHOT] [-m MAGNIFY]
                 [-g GSCRIPT] [-p PSCRIPT]
                 [-u POINTS_FNAME...]
                 [-q]
    pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]
                 [-s SHOT] [-m MAGNIFY]
                 [-g GSCRIPT] [-p PSCRIPT]
                 [-u POINTS_FNAME...]
                 [-q]
    pyNastranGUI -h | --help
    pyNastranGUI -v | --version

  Options:
    -h, --help                  show this help message and exit
    -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair,
                                             plot3d, stl, tetgen, usm3d)
    -i INPUT, --input INPUT     path to input file
    -o OUTPUT, --output OUTPUT  path to output file
    -g GSCRIPT, --geomscript GSCRIPT  path to geometry script file (runs before load geometry)
    -p PSCRIPT, --postscript PSCRIPT  path to post script file (runs after load geometry)
    -s SHOT, --shots SHOT       path to screenshot (only 1 for now)
    -m MAGNIFY, --magnify       how much should the resolution on a picture be magnified [default: 5]
    -u POINTS_FNAME, --user_points POINTS_FNAME               add user specified points to an alternate grid (repeatable)
    -q, --quiet                 prints debug messages (default=True)
    -v, --version               show program's version number and exit


The standard way to run the code:

.. code-block:: console

  >>> pyNastranGUI -f nastran -i model.bdf -o model1.op2 -o model2.op2

The **solid_bending.bdf** and **solid_bending.op2** files have been included
as examples that work in the GUI.  They are inside the "models" folder
(at the same level as setup.py).

Features
--------
 * Fringe Plot support
   - Custom Fringes supported
 * Command line interface
 * Scripting capability
 * High resolution screenshot (menu/button/keyboard)
 * Snap to Axis
 * Change Background Color

New Features
------------
 * nodal/element-based results at the same time
 * results may be shown alongside geometry
 * show/hide elements
   * can edit properties (e.g. color/opacity/size) using
     ``Edit Geometry Properties...`` on the ``View`` menu
   * additional points may be added with the ``-u`` option
 * attach custom CSV (comma-delimited) or .txt (space/tab-delimited) files as
   either node-based or element-based results
 * attach simplistic custom geometry
 * legend is more robust
 * clipping customization menu
 * save view menu
 * edges flippable from menu now

New Nastran Specific Features
-----------------------------
 * attach multiple OP2 files
 * supports SPOINTs
 * displacement/eigenvectors now shown as a deformation (real)
   * scale editable from legend menu
 * Edit Geometry Properties
   * SPC/MPC/RBE constraints
   * CAERO panel, subpanels
   * AEFACT control surfaces
   * SPLINE panels/points
   * bar/beam orientation vectors
   * CONM2

Supported Elements
------------------
 * CQUAD4 / CQUAD8
 * CTRIA3 / CTRIA6
 * CTETRA4 / CTETRA10
 * CHEXA8 / CHEXA20
 * CPENTA6 / CPENTA15
 * CSHEAR
 * CQUADR / CTRIAR
 * CBAR / CBEAM / CROD / CONROD / CELASx (displayed as lines)

BDF Requirements
----------------
 * Entire model can be cross-referenced
 * Same requirements as BDF (include an executive/case control deck, define
   all cross-referenced cards, etc.)

Scripting
---------
GUI commands are logged to the window with their call signature.  Users may
then use a custom Python script to take many pictures, show the sub-caero
panels, etc.  A sample CAERO script that shows individual CAERO subpanels
(instead of just the outline of the CAERO panel) is provided with the download.

For example, a model with CAERO elements:

.. image:: ../../pyNastran/gui/images/caero.png

can use the script:

.. code-block:: python

  # this script no longer works...
  self.show_cid = [1, 110000]
  self.is_sub_panels = True

and is called using:

.. code-block:: console

  >>> pyNastranGUI -f nastran -i bwb.bdf --script script.py


which creates:

.. image:: ../../pyNastran/gui/caero_subpanels.png

Scripting may be used to call any function in the GUI class.  This includes
dynamically loading geometry/results, changing results, taking screenshots,
rotating the model, etc.  Most of these commands are written to the
``COMMAND`` output.

Versioning Note
^^^^^^^^^^^^^^^
The GUI download is typically newer than the latest release version.

Additional Formats
------------------
Some of the results include:

   * **Nastran** ASCII input (\*.bdf, \*.nas, \*.dat, \*.pch); binary output (\*.op2)
     * geometry
        * node ID
        * element ID
        * property ID
        * material ID
        * thickness
        * normal
        * shell offset
        * PBAR/PBEAM/PBARL/PBEAML type

     * results (real only)
         * stress, strain
         * displacement, eigenvector, temperature, SPC forces, MPC forces, load vector

   * **Cart3d** ASCII/binary input (\*.tri); ASCII output (\*.triq)

     * Node ID
     * Element ID
     * Region
     * Cp, p, U, V, W, E, rho, rhoU, rhoV, rhoW, rhoE, Mach
     * Normal

   * **LaWGS** input (\*.wgs)

   * **Panair** input (\*.inp); output (agps, \*.out)

     * Patch ID
     * Normal X/Y/Z
     * Centroid X/Y/Z
     * Area
     * Node X/Y/Z
     * Cp

   * **STL** ASCII/binary input (\*.stl)

     * Normal X/Y/Z

   * **Tetgen** input (\*.smesh)

   * **Usm3d** surface input (\*.front, \*.cogsg); volume input (\*.cogsg); volume output (\*.flo)

     * Boundary Condition Region
     * Node ID
     * Cp, Mach, T, U, V, W, p, rhoU

Graphical Issues
^^^^^^^^^^^^^^^^
You'll have the best performance if you run the GUI on Windows with an new
NVIDIA graphics card and on a desktop.

If you're having issues, you should update the driver for your graphics card,
especially if you have a laptop or Radeon card. For a desktop machine, go to
the web site of the manufacturer of the graphics card. For a laptop, you should
normally go to the web site of the laptop manufacturer, though for NVIDIA you
may now find a newer driver available from NVIDIA.

Issues include:
  1. the backfaces of elements not being colored
  2. the GUI not working

Application Log Dislay
^^^^^^^^^^^^^^^^^^^^^^
It's possible that the Application Log will not be visible.  This is a PyQt4
install issue.  Reinstalling may fix the problem, but using one of the
recommended distributions is both the easiest and most reliable way to fix
this problem.



Features Overview
=================

Edit Geometry Properties
------------------------
The View -> "Edit Geometry Properties" menu brings up:

.. image:: ../../pyNastran/gui/images/edit_geometry_properties.png

This menu allows you to edit the opacity, line width, point size, show/hide various
things associated with the model.  The geometry auto-updates when changes are made.


Modify Legend
-------------
The View -> "Modify Legend" menu brings up:

.. image:: ../../pyNastran/gui/images/legend.png

This menu allows you to edit the max/min values of the legend as well as the orientation,
number format (e.g. float precision) and deflection scale.  Defaults are stored, so
they may always be gone back to.  The geometry will update when Apply/OK is clicked.
OK/Cancel will close the window.


Picking Results
---------------
Hover over an element and press the ``p`` key.  A label will appear.  This label will
appear at the centroid of an elemental result or the closest node to the selected location.
The value for the current result quantity will appear on the model.

.. image:: ../../pyNastran/gui/images/picking_results.png

For "NodeID", the xyz of the selcted point and the node in global XYZ space will be shown.
Labels may be cleared from the ``View`` menu.
Text color may also be changed from the ``View`` menu.


Focal Point
-----------
Hover over an element and press the ``f`` key.  The model will now rotate around
that point.


Model Clipping
--------------
Clipping let's you see "into" the model.

.. image:: ../../pyNastran/gui/images/clipping.png

Zoom in and hover over an element and press the ``f`` key.
The model will pan and now rotate around that point.
Continue to hold ``f`` while the model recenters.
Eventually, the frame will clip.
Reset the view by clicking the Undo-looking arrow.

Modify Groups
-------------
The View -> "Modify Groups" menu brings up:

.. image:: ../../pyNastran/gui/images/modify_groups1.png

Had you first clicked View -> "Create Groups by Property ID", you'd get:

.. image:: ../../pyNastran/gui/images/modify_groups2.png

Add/Remove use the "Patran-style" syntax:

.. code-block :: console

    # elements 1 to 10 inclusive
    1:10

    # elements 100 to the end
    100:#

    # every other element 1 to 11 - 1, 3, 5, 7, 9, 11
    # 1:11:2

The name of the group may also be changed, but duplicate names are not allowed.
The "main" group is the entire geometry.

The bolded/italicized text indicates the group that will be displayed to the screen.
The defaults will be updated when you click "Set As Main".  This will also update
the bolded/italicided group.

Camera Views
------------
The eyeball icon brings up a camera view.  You can set and save multiple camera views.
Additionally, views are written out for scripting.
You can script an external optimization process and take pictures every so many steps.

.. image:: ../../pyNastran/gui/images/camera_views.png


User Points
-----------

User points allow you to load a CSV of xyz points.
These may be loaded from within the GUI or from the command line.

.. code-block:: console

    # x, y, z
    1.0, 2.0, 3.0
    4.0, 5.0, 6.0

These will show up as points in the GUI with your requested filename.


User Geometry
-------------

User geometry is an attempt at creating a simple file format for defining geometry.
This may be loaded from the command line.  The structure will probably change.

.. code-block:: console

    # all supported cards
    #  - GRID
    #  - BAR
    #  - TRI
    #  - QUAD
    #
    # doesn't support:
    #  - solid elements
    #  - element properties
    #  - custom colors
    #  - coordinate systems
    #  - materials
    #  - loads
    #  - results

    #    id  x    y    z
    GRID, 1, 0.2, 0.3, 0.3
    GRID, 2, 1.2, 0.3, 0.3
    GRID, 3, 2.2, 0.3, 0.3
    GRID, 4, 5.2, 0.3, 0.3
    grid, 5, 5.2, 1.3, 2.3  # case insensitive

    #    ID, nodes
    BAR,  1, 1, 2
    TRI,  2, 1, 2, 3
    # this is a comment

    QUAD, 3, 1, 5, 3, 4
    QUAD, 4, 1, 2, 3, 4  # this is after a blank line


Custom Scalar Results
---------------------
Custom Elemental/Nodal CSV/TXT file results may be loaded.  The order and length is
important.  Results must be in nodal/elemental sorted order.  The following example
has 3 scalar values with 2 locations.

.. code-block:: console

      # x(%f), y(%i), z(%f)
      1.0,     2,     3.0
      4.0,     5,     6.0


Scripting
=========
The scripting menu allows for custom code and experimentation to be written without
loading a script from a file.  All valid Python is accepted.  Scripting commands should
start with ``self.``, but it's very powerful.  For example, you can do:

 - custom animations of mode shapes
 - high resolution screenshots
 - model introspection

High Resolution Screenshots
---------------------------

Option #1
^^^^^^^^^

.. code-block:: python

    self.on_take_screenshot('solid_bending.png', magnification=5)

Option #2
^^^^^^^^^

.. code-block:: python

    self.magnify = 5

Now take a screenshot.

Animation of Mode Shapes
------------------------

.. image:: ../../pyNastran/gui/images/solid_bending.gif

Attempt #1 - broken
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import time
    scales = [-1, 0.5, 0., 0.5, 1.0]
    title = 'Eigenvector'
    min_value = -1.0
    max_value = 1.0
    is_shown = True
    is_blue_to_red = True
    is_horizontal = False


    out = self.get_result_data_from_icase(icase)
    obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out

    # obj is NastranDisplacementResults
    min_value, max_value = obj.get_min_max(i, res_name)
    subtitle, label = self.get_subtitle_label(subcase_id)

    for scale in scales:
        #self.on_update_legend(title=title, min_value=min_value, max_value=max_value,
                              #scale=scale_value, data_format=data_format,
                              #is_blue_to_red=is_blue_to_red,
                              #is_discrete=is_discrete, is_horizontal=is_horizontal,
                              #is_shown=is_shown)
        self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
                                vector_size, subcase_id, result_type, location, subtitle, label,
                                revert_displaced=False)
        time.sleep(2)
        #self.grid.Modified()


Attempt #2 - broken
^^^^^^^^^^^^^^^^^^^


.. code-block:: python

    icase = 9
    out = self.get_result_data_from_icase(icase)
    obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out
    print(obj)

    scale = 1000.
    xyz = obj.xyz + scale * obj.dxyz[i, :]


Attempt #3 - broken
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    icase = 9
    out = self.get_result_data_from_icase(icase)
    obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out
    print(obj)

    label = ''
    min_value, max_value = obj.get_min_max(i, res_name)
    norm_value = float(max_value - min_value)
    scale = 100.

    name_vector = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
    case = obj
    xyz = obj.xyz + scale * obj.dxyz[i, :]

    grid_result_vector = self.set_grid_values(name_vector, case, vector_size, min_value, max_value, norm_value)


Attempt #3 - works
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from PIL.Image import open as open_image
    from pyNastran.gui.images2gif import writeGif

    icase = 9
    out = self.get_result_data_from_icase(icase)
    obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out

    xyz_base = obj.xyz
    nnodes = xyz_base.shape[0]
    actor = self.geometry_actors['main']

    screenshot_filenames = []
    scales = np.arange(-1., 1., 0.1) * 100.
    for scale in scales:
        screenshot_filename = 'solid_bending_%.0f.png' % scale
        xyz = xyz_base + scale * obj.dxyz[i, :]
        for j in range(nnodes):
            self.grid.GetPoints().SetPoint(j, xyz[j, :])

        self.grid.Modified()
        actor.Modified()
        self.rend.Render()
        self.on_take_screenshot(screenshot_filename, magnification=1)
        screenshot_filenames.append(screenshot_filename)
        break

    screenshot_filenames += screenshot_filenames[::-1][1:]
    gif_filename = 'solid_bending.gif'
    with open_image(screenshot_filenames[0]) as image:
        shape = (image.width, image.height)

    print('Writing gif to %s' % (gif_filename))
    shape2 = (shape[0] // 2, shape[1] // 2)
    images = [open_image(filename).resize(shape2) for filename in screenshot_filenames]  # this is wrong...

    #writeGif('solid_bending.gif', images, duration=1/framerate, subRectangles=False)
    writeGif(gif_filename, images, duration=0.1, dither=0)

