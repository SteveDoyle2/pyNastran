==============================
Graphical User Interface (GUI)
==============================

********
Overview
********

Introduction
============

The Graphical User Interface (GUI) looks like:

.. image:: ../../../pyNastran/gui/images/qt.png

A somewhat messy, but more featured image:

.. image:: ../../../pyNastran/gui/images/eigenvectors_groups_legend.png

The GUI also has a sidebar and transient support.


Purpose
=======
The pyNastran GUI was originally developed to solve a data validation problem.
It's hard to validate that things like coordinate systems were correct if you
can't look at the geometry in a native format.  As time went on, niche features
that were needed (e.g., aero panels) that were not supported natively in
Patran 2005, were added.  The goal is not to replace a code like Patran, but
complement it.

Since the intial development, the GUI has become significantly more capable
by adding features such as displacements and forces, so the need for a code like
Patran has decreased, but will not be eliminated.


Advantages of pyNastranGUI
--------------------------
 - command line interface for loading models
 - simple scripting
 - nice looking models
 - intuitive rotation
 - niche features
   - aero panels
   - aero splines
   - aero spline points
   - control surfaces
 - custom results from a CSV file
 - 64 bit support
   - Patran 2005 can't read in models that pyNastranGUI can
   - not an advantage for newer versions


Advantages of Patran
--------------------
 - CAD geometry support (e.g., IGES, Parasolid)
 - geometry creation (e.g., points, surfaces)
 - meshing
 - edit materials/properties
 - much better picking support
 - much better groups
 - better use of memory
 - grid point forces
 - many more...


Purpose of additional formats
=============================
Over time, pyNastran has also added converter and GUI support for additional
formats.  Nastran is not just once piece of the analysis puzzle and there is
a need for niche engineering formats.

While, you could convert a Cart3d model (a simple triangulation) to another
format like Nastran, you would need to map the geometry/result quantity of
interest (e.g., Mach Number) to something like Pressure.  That's unintuitive
and also requires writing an ill-defined format converter.  It's nice to load
it natively as you can also automatically create other quantities (e.g., the
bounding CFD box, free edges).

Finally, adding support for alternate formats drives GUI development.  The model
reload functionality was added to address loading the latest time step of a Usm3d
model.  It was repurposed to reload the geometry for other formats.  This is very
useful when creating aero panels and you want to see your changes.  The groups
functionality benefits all formats.


Additonal formats include:

   - panair
   - cart3d
   - stl
   - tecplot
   - AFLR
    - bsurf
    - surf
    - ugrid
   - usm3d


Setup Note
==========
Download the entire package from Github or just the `GUI
<https://sourceforge.net/projects/pynastran/files/?source=navbar/>`_ executable.

If you download the source, make sure you follow the `Installation Guide
<https://github.com/SteveDoyle2/pyNastran/wiki/Installation>`_ and use
**setup.py develop** and not **setup.py install**.

Python 2.7 with ``vtk==5.10.1`` or ``vtk==6.3.0`` will give you the best looking
GUI.

VTK 7 issues:
 - edges can't be colored black.

Python 3 issues:
 - the GUI in Python 3 won't save your settings (e.g., background color).


Running the GUI
===============
On the command line:

.. code-block:: console

  >>> pyNastranGUI

To view the options:

.. code-block:: console

  >>> pyNastranGUI --help

    Usage:
      pyNastranGUI [-f FORMAT] INPUT [-o OUTPUT]
                   [-s SHOT] [-m MAGNIFY]
                   [-g GSCRIPT] [-p PSCRIPT]
                   [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]
                   [-q] [--groups]
      pyNastranGUI [-f FORMAT] INPUT OUTPUT [-o OUTPUT]
                   [-s SHOT] [-m MAGNIFY]
                   [-g GSCRIPT] [-p PSCRIPT]
                   [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]
                   [-q] [--groups]
      pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]
                   [-s SHOT] [-m MAGNIFY]
                   [-g GSCRIPT] [-p PSCRIPT]
                   [-u POINTS_FNAME...] [--user_geom GEOM_FNAME...]
                   [-q] [--groups]
      pyNastranGUI -h | --help
      pyNastranGUI -v | --version

    Primary Options:
      -f FORMAT, --format FORMAT  format type (avus, cart3d, lawgs, nastran, panair,
                                               plot3d, stl, surf, tetgen, ugrid, usm3d)
      -i INPUT, --input INPUT     path to input file
      -o OUTPUT, --output OUTPUT  path to output file

    Secondary Options:
      -g GSCRIPT, --geomscript        path to geometry script file (runs before load geometry)
      -p PSCRIPT, --postscript        path to post script file (runs after load geometry)
      -s SHOT, --shots SHOT           path to screenshot (only 1 for now)
      -m MAGNIFY, --magnify           how much should the resolution on a picture be magnified [default: 5]
      --groups                        enables groups
      --user_geom GEOM_FNAME          add user specified points to an alternate grid (repeatable)
      -u POINTS_FNAME, --user_points  add user specified points to an alternate grid (repeatable)

    Info:
      -q, --quiet    prints debug messages (default=True)
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit



The standard way to run the code:

.. code-block:: console

  >>> pyNastranGUI -f nastran -i model.bdf -o model1.op2 -o model2.op2

The **solid_bending.bdf** and **solid_bending.op2** files have been included
as examples that work in the GUI.  They are inside the "models" folder
(at the same level as ``setup.py``).

You can also run it like:

  >>> pyNastranGUI model.bdf model1.op2

Here the code will guess based on your file extension what your file format is.
If you want to load a second OP2, you must use ``-o model2.op2``.

Features
========
 * fringe plot support

   * elemental/nodal results
   * custom CSV results

 * deflection results
 * force results

 * command line interface
 * scripting capability
 * high resolution screenshot
 * show/hide elements

   * can edit properties (e.g. color/opacity/size) using
     ``Edit Geometry Properties...`` on the ``View`` menu

 * legend menu
 * save/load view menu

Minor Features
==============
 * snap to axis
 * clipping customization menu
 * edges flippable from menu
 * change label color/size menu
 * change background color
 * attach simplistic custom geometry
 * additional points may be added with the ``-u`` option


Nastran Specific Features
=========================
 * attach multiple OP2 files
 * supports SPOINTs
 * displacement/eigenvectors now shown as a deformation (real)

   * scale editable from legend menu
   * rotated into global frame

 * Edit Geometry Properties

   * SPC/MPC/RBE constraints
   * CAERO panel, subpanels
   * AEFACT control surfaces
   * SPLINE panels/points
   * bar/beam orientation vectors
   * CONM2


BDF Requirements
================
 * Entire model can be cross-referenced
 * Same requirements as BDF (include an executive/case control deck, define
   all cross-referenced cards, etc.)


Versioning Note
---------------
The GUI download is typically newer than the latest release version.

Additional Formats
==================
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
        * element quality (min/max interior angle, skew angle, taper ratio, area ratio)

     * real results
         * stress, strain
         * displacement, eigenvector, temperature, SPC forces, MPC forces, load vector
     * complex results
         * displacement, eigenvector

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


*****************
Features Overview
*****************

Edit Geometry Properties
========================
The View -> "Edit Geometry Properties" menu brings up:

.. image:: ../../../pyNastran/gui/images/edit_geometry_properties.png

This menu allows you to edit the opacity, line width, point size, show/hide various
things associated with the model.  The geometry auto-updates when changes are made.


Modify Legend
=============
The View -> "Modify Legend" menu brings up:

.. image:: ../../../pyNastran/gui/images/legend.png

This menu allows you to edit the max/min values of the legend as well as the orientation,
number format (e.g. float precision) and deflection scale.  Defaults are stored, so
they may always be gone back to.  The geometry will update when Apply/OK is clicked.
OK/Cancel will close the window.

Animation of Displacment/Mode Shapes
====================================

The animation menu is a sub-menu found on the Legend Menu.  It supports:
 - Scale Factor
 - Total Time (sec)
 - Frames/Second
 - Resolution Scale
 - Output Directory
 - Filename
 - scale/phase animation
 - one/two sided animations
 - infinite looping
 - delete images
 - make gif
 - TODO: animate time/frequency/loadstep
 - TODO: progate results based on result locking

You must load the animation menu when a displacement-like result is active.
You may then change to a scalar result to show during the animation.


Picking Results
===============
Hover over an element and press the ``p`` key.  A label will appear.  This label
will appear at the centroid of an elemental result or the closest node to the
selected location.  The value for the current result quantity will appear on the
model.  You may also use the button.

.. image:: ../../../pyNastran/gui/images/picking_results.png

For "NodeID", the xyz of the selcted point and the node in global XYZ space will be shown.
Labels may be cleared from the ``View`` menu.
Text color may also be changed from the ``View`` menu.


Focal Point
===========
Hover over an element and press the ``f`` key.  The model will now rotate around
that point.  You may also use the button.


Model Clipping
==============
Clipping let's you see "into" the model.

.. image:: ../../../pyNastran/gui/images/clipping.png

Zoom in and hover over an element and press the ``f`` key.
The model will pan and now rotate around that point.
Continue to hold ``f`` while the model recenters.
Eventually, the frame will clip.
Reset the view by clicking the Undo-looking arrow at the top.

Modify Groups
=============
The View -> "Modify Groups" menu brings up:

.. image:: ../../../pyNastran/gui/images/modify_groups1.png

Had you first clicked View -> "Create Groups by Property ID", you'd get:

.. image:: ../../../pyNastran/gui/images/modify_groups2.png

Add/Remove use the "Patran-style" syntax:

.. code-block :: console

    # elements 1 to 10 inclusive
    1:10

    # elements 100 to the end
    100:#

    # every other element 1 to 11 - 1, 3, 5, 7, 9, 11
    1:11:2

The name of the group may also be changed, but duplicate names are not allowed.
The "main" group is the entire geometry.

The bolded/italicized text indicates the group that will be displayed to the screen.
The defaults will be updated when you click ``Set As Main``.  This will also update
the bolded/italicided group.

Camera Views
============
The eyeball icon brings up a camera view.  You can set and save multiple camera views.
Additionally, views are written out for scripting.
You can script an external optimization process and take pictures every so many steps.

.. image:: ../../../pyNastran/gui/images/camera_views.png


User Points
===========

User points allow you to load a CSV of xyz points.
These may be loaded from within the GUI or from the command line.

.. code-block:: console

    # x, y, z
    1.0, 2.0, 3.0
    4.0, 5.0, 6.0

These will show up as points in the GUI with your requested filename.


User Geometry
=============

User geometry is an attempt at creating a simple file format for defining geometry.
This may be loaded from the command line.  The structure will probably change.

The geometry may be modified from the ``Edit Geometry Properties`` menu.

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
=====================
Custom Elemental/Nodal CSV/TXT file results may be loaded.  The order and
length is important.  Results must be in nodal/elemental sorted order.
The following example has 3 scalar values with 2 locations.  The model must
have **only** two nodes.

.. code-block:: console

      # x(%f), y(%i), z(%f)
      1.0,     2,     3.0
      4.0,     5,     6.0

Custom Results Specific Buttons
===============================
Nastran Static/Dynamic Aero solutions require custom cards that create
difficult to view, difficult to validate geometry.  The pyNastranGUI
aides in creating models.  The CAERO panels are seen when a model is loaded:

.. image:: ../../../pyNastran/gui/images/caero.png

Additionally, by clicking the ``Toggle CAERO Subpanels`` button
(the figure is somewhat outdated), the subpanels may be seen:

.. image:: ../../../pyNastran/gui/images/caero_subpanels.png

Additionally, flaps are shown from within the GUI.  SPLINE surfaces
are also generated and may be seen on the ``View`` -> ``Edit Geometry Properties``
menu.

*********
Scripting
*********
GUI commands are logged to the window with their call signature.
Scripting may be used to call any function in the GUI class.
Most of these commands are written to the ``COMMAND`` output.

For example, you can:

 - load geometry
 - load results
 - plot unsupported result types
 - custom animations of mode shapes
 - high resolution screenshots
 - model introspection


Using the scripting menu
========================
The scripting menu allows for custom code and experimentation to be written
without loading a script from a file.  All valid Python is accepted.
Scripting commands should start with ``self.`` as they're left off from the menu.
Local variables do not need this.


Command line scripting
======================
``geom_script`` runs after the load_geometry method, while
``postscript`` runs after load_results has been performed

.. code-block :: python

    import sys
    self.on_take_screenshot('solid_bending.png', magnify=5)
    sys.exit()

.. code-block :: console

    >>> pyNastranGUI solid_bending.bdf solid_bending.op2 --postscript take_picture.py


High Resolution Screenshots
===========================

Option #1
---------

.. code-block:: python

    self.on_take_screenshot('solid_bending.png', magnify=5)

Option #2
---------

.. code-block:: python

    self.magnify = 5

Now take a screenshot.

Animation of Displacment/Mode Shapes
====================================

While it's possible to take multiple screenshots of geometry with
different scale factors, it's tedious.  Additionally, you can only
plot displacement-type results (e.g. displacement, eigenvector)
with deflection and not result types like Node ID or stress
unless you write a script.

.. image:: ../../../pyNastran/gui/images/solid_bending.gif


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
        self.on_take_screenshot(screenshot_filename, magnify=1)
        screenshot_filenames.append(screenshot_filename)

    screenshot_filenames += screenshot_filenames[::-1][1:]
    gif_filename = 'solid_bending.gif'
    with open_image(screenshot_filenames[0]) as image:
        shape = (image.width, image.height)

    print('Writing gif to %s' % (gif_filename))

    # down-res the image so we use less space
    shape2 = (shape[0] // 2, shape[1] // 2)
    images = [open_image(filename).resize(shape2) for filename in screenshot_filenames]

    #writeGif('solid_bending.gif', images, duration=1/framerate, subRectangles=False)
    writeGif(gif_filename, images, duration=0.1, dither=False)

.. Attempt #2 - broken
.. -------------------
..
.. .. code-block:: python
..
..     import time
..     scales = [-1, 0.5, 0., 0.5, 1.0]
..     title = 'Eigenvector'
..     min_value = -1.0
..     max_value = 1.0
..     is_shown = True
..     is_blue_to_red = True
..     is_horizontal = False
..
..
..     out = self.get_result_data_from_icase(icase)
..     obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out
..
..     # obj is NastranDisplacementResults
..     min_value, max_value = obj.get_min_max(i, res_name)
..     subtitle, label = self.get_subtitle_label(subcase_id)
..
..     for scale in scales:
..         #self.on_update_legend(title=title, min_value=min_value, max_value=max_value,
..                               #scale=scale_value, data_format=data_format,
..                               #is_blue_to_red=is_blue_to_red,
..                               #is_discrete=is_discrete, is_horizontal=is_horizontal,
..                               #is_shown=is_shown)
..         self._final_grid_update(name_vector, grid_result_vector, obj, i, res_name,
..                                 vector_size, subcase_id, result_type, location, subtitle, label,
..                                 revert_displaced=False)
..         time.sleep(2)
..         #self.grid.Modified()


.. Attempt #3 - broken
.. -------------------
..
.. .. code-block:: python
..
..     icase = 9
..     out = self.get_result_data_from_icase(icase)
..     obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out
..     print(obj)
..
..     label = ''
..     min_value, max_value = obj.get_min_max(i, res_name)
..     norm_value = float(max_value - min_value)
..     scale = 100.
..
..     name_vector = (vector_size, subcase_id, result_type, label, min_value, max_value, scale)
..     case = obj
..     xyz = obj.xyz + scale * obj.dxyz[i, :]
..
..     grid_result_vector = self.set_grid_values(name_vector, case, vector_size, min_value, max_value, norm_value)


.. Animation of Complex Mode Shapes
.. ================================


.. Complex Mode Shapes (not done)
.. ------------------------------

 .. code-block:: python

     from PIL.Image import open as open_image
     from pyNastran.gui.images2gif import writeGif

     from pyNastran.op2.op2 import read_op2
     model = read_op2(op2_filename)

     xyz_undef = self.xyz_cid0
     nnodes = xyz_undef.shape[0]

     #out = self.get_result_data_from_icase(icase)
     #obj, i, j, res_name, subcase_id, result_type, vector_size, location, data_format, label2 = out
     actor = self.geometry_actors['main']

     subcase_id = 1
     imode = 10
     eigenvectors = model.eigenvectors[subcase_id].data[imode - 1,:,:]

     #-------------------------------------------------------------------
     mag = np.abs(eigenvectors[:, :3])
     phase = np.angle(eigenvectors[:, :3])
     reals = np.real(eigenvectors[:, :3])
     imags = np.imag(eigenvectors[:, :3])

     nframes = 10
     amplitude = np.ones(nframes) * 5
     screenshot_filenames = []
     for i in range(nframes):
         screenshot_filename = 'solid_bending_complex_%i.png' % i

         theta = (2*np.pi * i/nframes) % (2*np.pi)
         defl = amplitude[i] * (reals*np.cos(theta) + imags*np.sin(theta))
         xyz_def = xyz_undef + defl
         for j in range(nnodes):
             self.grid.GetPoints().SetPoint(j, xyz_def[j, :])

         self.grid.Modified()
         actor.Modified()
         self.rend.Render()
         self.on_take_screenshot(screenshot_filename, magnify=1)
         screenshot_filenames.append(screenshot_filename)
     screenshot_filenames += screenshot_filenames[::-1][1:]


     #-------------------------------------------------------------------

     gif_filename = 'solid_bending_complex.gif'
     with open_image(screenshot_filenames[0]) as image:
         shape = (image.width, image.height)

     print('Writing gif to %s' % (gif_filename))

     # down-res the image so we use less space
     shape2 = (shape[0] // 2, shape[1] // 2)
     images = [open_image(filename).resize(shape2) for filename in screenshot_filenames]

     writeGif(gif_filename, images, duration=0.1, dither=False)
