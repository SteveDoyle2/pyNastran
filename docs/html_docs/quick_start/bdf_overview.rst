============
BDF Overview
============

------------
Introduction
------------
This is meant as a tutorial on how to use the pyNastran ``pyNastran.bdf.bdf.BDF`` class


The **head**/**tail**/**file_slice** methods can be found at:

    https://github.com/SteveDoyle2/pyNastran/blob/v0.7/docs_sphinx/manual/py_docs/bdf_doc.py

These examples can be found at:

    https://github.com/SteveDoyle2/pyNastran/blob/v0.7/docs_sphinx/manual/py_docs/bdf_doc.py

---------------------
Example 1: Read/Write
---------------------
this example will demonstate:

 - reading the BDF

 - getting some basic information

 - writing the BDF

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> bdf_filename = os.path.join(test_path, 'solid_bending.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename)

For unicode:

  The standard encoding is utf-8, but most English decks should use latin1 and will fail with utf-8.

  If you just have ascii, then you don't need to worry about the encoding.

.. code-block:: python

    >>> model.read_bdf(bdf_filename, encoding='latin1')

print information about the model

.. code-block:: python

    >>> print(model.get_bdf_stats())
    ---BDF Statistics---
    SOL 101

    bdf.loads[1]
      FORCE:   23

    bdf.loads[2]
      LOAD:    1

    bdf.params
      PARAM    : 2

    bdf.nodes
      GRID     : 72

    bdf.elements
      CTETRA   : 186

    bdf.properties
      PSOLID   : 1

    bdf.materials
      MAT1     : 1

    bdf.coords
      CORD2R   : ???

write the file

.. code-block:: python

    >>> bdf_filename_out = os.path.join(test_path, 'solid_bending_out.bdf')
    >>> model.write_bdf(bdf_filename_out)

looking at the output

.. code-block:: python

    >>> print(file_slice(bdf_filename_out, 94, 100))
    GRID          71         .500008 1.61116      3.
    GRID          72         .500015 1.00001      3.
    $ELEMENTS_WITH_PROPERTIES
    PSOLID         1       1
    CTETRA         1       1       8      13      67      33
    CTETRA         2       1       8       7      62      59

write the file with large field format; double precision

.. code-block:: python

    >>> bdf_filename_out2 = os.path.join(test_path, 'solid_bending_out2.bdf')
    >>> model.write_bdf(bdf_filename_out2, size=16, is_double=False)
    >>> print(file_slice(bdf_filename_out2, 166, 175))
    GRID*                 71                         .500008         1.61116
    *                     3.
    GRID*                 72                         .500015         1.00001
    *                     3.
    $ELEMENTS_WITH_PROPERTIES
    PSOLID         1       1
    CTETRA         1       1       8      13      67      33
    CTETRA         2       1       8       7      62      59
    CTETRA         3       1       8      45      58      66

write the file with large field format; double precision

.. code-block:: python

    >>> bdf_filename_out3 = os.path.join(test_path, 'solid_bending_out3.bdf')
    >>> model.write_bdf(bdf_filename_out3, size=16, is_double=True)
    >>> print(file_slice(bdf_filename_out3, 166, 175))
    GRID*                 71                5.0000800000D-011.6111600000D+00
    *       3.0000000000D+00
    GRID*                 72                5.0001500000D-011.0000100000D+00
    *       3.0000000000D+00
    $ELEMENTS_WITH_PROPERTIES
    PSOLID         1       1
    CTETRA         1       1       8      13      67      33
    CTETRA         2       1       8       7      62      59
    CTETRA         3       1       8      45      58      66

--------------------------
Example 2:  Printing Nodes
--------------------------
this example will demonstate:

 - writing cards

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> bdf_filename = os.path.join(test_path, 'solid_bending.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

Method 1 - using objects
------------------------


GRIDs

.. code-block:: python

    >>> for nid,node in sorted(model.nodes.items()):
    >>>     f.write(node.write_card(size=8, is_double=False))

GRIDSET

.. code-block:: python

    >>> if model.gridSet:
    >>>     f.write(model.gridSet.write_card(size=8, is_double=False))

SPOINTs

.. code-block:: python

    >>> if model.spoints:
    >>>     f.write(model.spoints.write_card(size=8, is_double=False))

Coord

.. code-block:: python

    >>> for cid,coord in sorted(model.coords.items()):
    >>>     if cid != 0:  # if CID=0 is the global frame, skip it
    >>>         f.write(coord)

Method 2 - using built-in methods
---------------------------------

    >>> model._write_nodes(f)
    >>> model._write_coords(f)

----------------------------------------
Example 3:  Printing Elements/Properties
----------------------------------------
Print the Element ID and associated Node and Property to an Output File

note this skips rigidElements

this example will demonstate:

 - using the BDF class to write cards/properties

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> bdf_filename = os.path.join(test_path, 'solid_bending.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

Method 1 - using objects
------------------------

    >>> for eid, element in sorted(model.elements.items()):
    >>>     f.write(element.write_card(size=8, is_double=False))
    >>> for pid, prop in sorted(model.properties.items()):
    >>>     f.write(prop.write_card(size=8, is_double=False))

Method 2 - using built-in method
--------------------------------

    >>> model._write_elements_properties(f)

Method 3 - using built-in methods
---------------------------------

    >>> model._write_elements(f)
    >>> model._write_properties(f)

--------------------------------
Example 4: Get Element ID & Type
--------------------------------
Print the Element ID and its type(e.g. CQUAD4, CTRIA3, etc.) to a file

note this skips rigidElements

this example will demonstate:

 - accessing element type information

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> bdf_filename = os.path.join(test_path, 'solid_bending.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

Method 1 - using objects
---------------------------------

    >>> for eid,element in sorted(model.elements.items()):
    >>>     msg = 'eid=%s type=%s\n' %(eid, element.type)
    >>> f.write(msg)

----------------------------------
Example 5: Get Elements by Node ID
----------------------------------
this example will demonstate:

 - getting the list of elements that share a certain node

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'solid_bending')
    >>> bdf_filename = os.path.join(test_path, 'solid_bending.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

given a Node, get the Elements Attached to that Node

assume node 55

doesnt support 0d/1d elements yet

.. code-block:: python

    >>> nid_to_eids_map = model.get_node_id_to_element_ids_map()
    >>> eids = nid_to_eids_map[55]

convert to elements instead of element IDs

.. code-block:: python

    >>> elements = []
    >>> for eid in eids:
    >>>     elements.append(model.Element(eid))
    >>> print("eids = %s" % eids)
    >>> print("elements =\n %s" % elements)

---------------------------------------
Example 6:  Get Elements by Property ID
---------------------------------------
this example will demonstate:

 - getting a list of elements that have a certain property

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'sol_101_elements')
    >>> bdf_filename = os.path.join(test_path, 'static_solid_shell_bar.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

Creating a List of Elements based on a Property ID

assume pid=1

.. code-block:: python

    >>> pid_to_eids_map = model.get_property_id_to_element_ids_map()
    >>> eids4  = pid_to_eids_map[4] # PSHELL
    >>> print("eids4 = %s" % eids4)
    eids4 = [6, 7, 8, 9, 10, 11]

convert to elements instead of element IDs

.. code-block:: python

    >>> elements4 = []
    >>> for eid in eids4:
    >>>     elements4.append(model.Element(eid))

just to verify

.. code-block:: python

    >>> elem = model.elements[eids4[0]]
    >>> print(elem.pid)
    PSHELL         4       1     .25       1               1

---------------------------------------
Example 7:  Get Elements by Material ID
---------------------------------------
this example will demonstate:

 - getting a list of elements that have a certain material

our model

.. code-block:: python

    >>> import pyNastran
    >>> pkg_path = pyNastran.__path__[0]
    >>> test_path = os.path.join(pkg_path, '..', 'models', 'sol_101_elements')
    >>> bdf_filename = os.path.join(test_path, 'static_solid_shell_bar.bdf')

instantiate the model

.. code-block:: python

    >>> from pyNastran.bdf.bdf import BDF
    >>> model = BDF()
    >>> model.read_bdf(bdf_filename, xref=True)
    >>> f = open('junk.out', 'w')

assume you want the eids for material 10

.. code-block:: python

    >>> pid_to_eids_map = model.get_property_id_to_element_ids_map()
    >>> mid_to_pids_map = model.get_material_id_to_property_ids_map()
    >>> pids1 = mid_to_pids_map[1]
    >>> print('pids1 = %s' % pids1)
    pids1 = [1, 2, 3, 4, 5]
    >>> eids = []
    >>> for pid in pids1:
    >>>     eids += pid_to_eids_map[pid]

convert to elements instead of element IDs

.. code-block:: python

    >>> elements = []
    >>> for eid in eids:
    >>>     element = model.Element(eid)
    >>>     elements.append(element)
    >>>     print(str(element).rstrip())

    CBAR          13       1      15      19      0.      1.      0.
    $ Direct Text Input for Bulk Data
    $ Pset: "shell" will be imported as: "pshell.1"
    CHEXA          1       2       2       3       4       1       8       5
                   6       7
    CPENTA         2       2       6       8       5      10      11       9
    CPENTA         3       2       6       7       8      10      12      11
    CTETRA         4       2      10      11       9      13
    CTETRA         5       2      10      12      11      13
    CROD          14       3      16      20
    CROD          15       3      17      21
    CQUAD4         6       4       4       1      14      15
    CQUAD4         7       4       3       2      17      16
    CTRIA3         8       4       4       3      16
    CTRIA3         9       4      16      15       4
    CTRIA3        10       4       1       2      17
    CTRIA3        11       4      17      14       1
    $
    CBEAM         12       5      14      18      0.      1.      0.     GGG
