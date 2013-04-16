    
==============================
Developer: BDF Reading
==============================

This document is a reference for developers of pyNastran, but are not necessarily
 people that are familiar with Nastran or even Finite Element Analysis (FEA).

:mod:`bdf`:   Introduction
----------------------------

:mod:`bdf` module controls the model object that is instantiated with  
``model = BDF()`` the ``BDF.__init__`` method
is called when ``model = BDF()`` is run.  The :py:class:`pyNastran.bdf.bdf.BDF` 
used to be very large and was split (in a non-standard way) into multiple files 
that allowed for simpler development.  The classes:

 * :py:class:`pyNastran.bdf.bdfInterface.bdf_Reader.BDFReader`,
 * :py:class:`pyNastran.bdf.bdfInterface.getCard.GetMethods`,
 * :py:class:`pyNastran.bdf.bdfInterface.addCard.AddMethods`, 
 * :py:class:`pyNastran.bdf.bdfInterface.bdf_writeMesh.WriteMesh`,
 * :py:class:`pyNastran.bdf.bdfInterface.bdf_cardMethods.CardMethods`,
 * :py:class:`pyNastran.bdf.bdfInterface.crossReference.XrefMesh` 
 
 
are basically bags of functions for the "model" object.

The :py:attr:`pyNastran.bdf.bdf.BDF.cardsToRead` attribute limits the cards that 
:mod:`pyNastran` processes and can be modified by the user in order to fix bugs 
or make their code run faster.
     
Moving onto :py:attr:`pyNastran.bdf.bdf.BDF._init_solution` sets a series of 
alternate names for Nastran solution types.  For example, a solution 101 is 
a static solution (no acceleration) and will calculate the displacements of 
the system :math:`[K]\{x\} = \{F\}`.  You can then extract stresses and strains.
Other solution numbers solve different equations.

In methods :py:meth:`pyNastran.bdf.bdf.BDF._init_structural_defaults`, 
:py:meth:`pyNastran.bdf.bdf.BDF._init_aero_defaults`,
:py:meth:`pyNastran.bdf.bdf.BDF._init_thermal_defaults` 
the structural, aerodyanmic, and thermal card holders (e.g. model.materials)
are defined as dictionaries.

Finally, the :py:meth:`pyNastran.bdf.bdf.BDF.readBDF` method is defined.  
There are three sections to a BDF/DAT/NAS file. BDF (Bulk Data File) is the file
format used by MSC Nastran and the DAT (data?) file is used by NX Nastran.
NAS (Nastran) is also used.

The first section is the "Executive Control Deck" and contains a "SOL 101" or 
"SOL 103" or "SOL STATIC" depending on the solution type. It ends when the "CEND"
marker is found. Then the "Case Control Deck" is read. Here, general solution
information is listed, such  as what outputs do you want and what loads are applied
(not all loads in the file are necessarily applied).  Finally this section defines
one or more subcases, which are different load combinations. The last section 
is the "Bulk Data Deck".  This section stores 99% of the file and this section 
introduces very specific formatting restrictions for "cards".

A basic model will be made of nodes, elements, properties, and materials.  
For example, for a square plate made of steel model GRID cards are used for
the nodes, CQUAD4 or CTRIA3 cards are used for the elements (the C generally 
indicates the card is an element so quadrilateral element and triangular element).
The element has a property (PSHELL) that defines the thickness.  Similarly, 
properties generally start with P.  Finally,  the property references a material
(MAT1) that defines the material as steel.  INCLUDE cards may also be used to 
add additional files into the BDF.


:mod:`bdf`: Card Formatting
-----------------------------

A "card" is at most 72-characters wide.  Comments may follow the card if a 
$ sign is used.

The standard card is called small field format (single precision) and has 9 fields
defined per line, each with 8-characters and are fixed width.  Multiline cards are 
implied by leaving 8 spaces at the beginning of the following line.  
Alternatively, a + sign may be used in the first 8 spaces.

The large field format (double precision) card uses a :math:`1 * 8 + 4 * 16`
to reach the 72-character width instead of :math:`1 * 8 + 8 * 8` characters.  
If the first line of a card is double precision, a * follows the card name, 
so all card  names are 7-characters or less.  If the second line of a card is 
double precision, a * begins the line.  A single line of a small field formatted
takes exactly two lines to write if large field format is used.

The CSV (comma separated value) format is similar to small field format.  
It's less picky then the 8-character format, but much harder to read.  
It is still subject to the 9 fields per line restriction.  If a CSV card has 
a * on it, the card becomes a large field CSV formatted card and may have only 
5 fields on the line (including the blank field).

Although introduced as separate types of cards, small field format and large
field format may be mixed and matched. However, this only occurs for hand-edited 
BDFs.  There's also a difficult to understand format known as a "continuation card".
This uses values from previous cards and is basically a *for* loop.
Hundreds of cards may be defined in just a few lines.



:mod:`bdf` : Parsing
----------------------

A basic card is a GRID card.  Once parsed, a standard grid card will have fields
of ``['GRID', nodeID, coord, x, y, z]``. This section will discuss how a card is
parsed.

The :py:meth:`pyNastran.bdf.bdf.BDF.readBDF` method must generalize the way
files are opened because INCLUDE files may be used. Once the Executive and Case
Control Decks are read (using simple while loops), the
:py:meth:`pyNastran.bdf.bdf.BDF._read_bulk_data_deck` method is called.

This method (:meth:`BDF._read_bulk_data_deck`) keeps looping over the file as
long as there are files opened (an INCLUDE file side effect) and calls:
``(rawCard, card, cardName) = self._get_card(debug=False)``. ``cardName`` is
just the card's name, while ``rawCard`` is the full, unparsed card. ``card`` is
the card  parsed into fields as a ``card`` object, which is basically a list
of fields ``['GRID', nodeID, coord, x, y, z]``.

The basic idea of the ``self._get_card()`` method(see
:py:meth:`pyNastran.bdf.bdfInterface.bdf_cardMethods.CardMethods._get_card`) is
to make a ``self.linesPack``, (:py:attr:`pyNastran.bdf.bdf.BDF.linesPack`)which
is a list of 1500 lines that are stored in memory.  When a new line is added to
``self.linesPack``, it is first stripped of comments in order to simplify
parsing.  If the linesPack is empty or 50 blank lines are found, the code
assumes an infinite loop has been entered and closes the file.  If additional
files are open, the ``linesPack`` from the previous file will be used (INCLUDE
file case).

Now that we have 1500 lines in linesPack, we must call:  ``(i, tempcard) =
self._get_multi_line_card(i, tempcard)`` to get the card. ``tempcard`` starts
out as the first line in the card and afterwards contains all lines of the card.
``tempcard`` will eventually become ``rawCard``.  It's important to note the
``tempcard`` passed into
:py:meth:`pyNastran.bdf.bdfInterface.bdf_cardMethods.CardMethods._get_multi_line_card`
is a 72-character string (generally) and the ``tempcard`` output is a list of
8-character (or 16) width fields.  Why: the original data isn't needed, so the
variable is reused.

:py:meth:`pyNastran.bdf.bdfInterface.bdf_cardMethods.CardMethods._get_multi_line_card`
will search through the ``linesPack`` and try to end the card by looking for a
non-whitespace character in the first character position (all cards must have
``field[0]`` justified).
If a * is found, it's double precision, if not it's small field.  Additionally,
if a ',' is found it's CSV format.
So the formats are:

 #. small field, 
 #. small field CSV, 
 #. large field, 
 #. large field CSV.


Once the format of the line is known, it's an easy process to split the card
(see
:py:meth:`pyNastran.bdf.bdfInterface.bdf_cardMethods.CardMethods.processCard`)
and turn it into a :py:class:`pyNastran.bdf.bdfInterface.BDF_Card.BDFCard`
object.  Note that the first field in any line beyond the first one must be
blank and is ignored.  This prevents cards defined in small field and large
field to have the same value defined in different positions of the list.

Finally, as Nastran is very specific in putting a decimal on float values, it's
easy to parse field values into their proper type dynamically.  This is
especially important when a field may be defined as an integer, float, a string,
or be left blank and the variable is different depending on variable type. 
Strings, must being with alphabetical characters (A, B, C) and are case
insensitive, which is why a "GRID" card is called a "GRID" card and not a "grid"
card.


:mod:`bdf` : Card Object
--------------------------

A :py:class:`pyNastran.bdf.bdfInterface.BDF_Card.BDFCard` object is basically a
list of fields of ``['GRID', nodeID, coord, x, y, z]`` with methods to get the
1st entry (``nodeID``) as ``card.field(1)`` instead of ``fields[1]`` for a list.
A card object is useful for setting defaults.  The ``x, y``, and ``z`` values
on the GRID card have defaults of 0.0, so ``card.field(3,0.0)`` may be used to
get the ``x`` coordinate. Finally, ``card.fields(3,5,[0.,0.,0.])`` may be used
to get ``xyz`` and set the defaults in a single step.  Additionally, the
``card`` object is useful when parsing "continuation cards", but is typically
disabled.

After an excessively long branch of ``cardNames`` in
:py:meth:`pyNastran.bdf.BDF.readBDF`, the card object is turned into a GRID,
CTRIA3, CQUAD4, PSHELL, MAT1 or any of 200 other card types.  There are roughly
as many nodes as there are elements, which make up roughly 95% of the cards in
large models.  The difference in a large model and a small model, is the
discretization and will change nodes, elements, loads, and constraints.  Loads
and constraints are applied to only small portions of the model and (generally)
only the boundary of a model.  The number of propertie and materials is very
likely the same.

Most cards are stored in a dictionary based on their integer ID.  IDs may be
used only once, but if a card is exactly duplicated, it is still valid.


:mod:`shell`: CQUAD4 Object
-----------------------------
In ``bdf/cards/elements/shell.py``, the
:py:class:`pyNastran.bdf.cards.elements.shell.CQUAD4` is defined.

The :py:class:`pyNastran.bdf.cards.elements.shell.CQUAD4` is a shell-type
element and must reference a PSHELL (isotropic property) or a PCOMP (composite
property) card.  An example of an isotropic material is steel or aluminum and a
composite material would be fiberglass or layers of carbon fiber/epoxy resin at
layed up at different angles.

The PSHELL may reference MAT1 (isotropic material) cards, while the PCOMP card
may reference MAT1 or MAT8 (orthotropic material) cards.  An orthotropic
material is stronger in longitudinally than laterally (e.g. fibers are oriented
unidirectionally in a carbon fiber composite).

The :py:class:`pyNastran.bdf.cards.elements.shell.CQUAD4` class inherits from
the :py:class:`pyNastran.bdf.cards.elements.shell.QuadShell` class which defines
common methods to the various QUAD-type cards.  There are additional QUAD
element with different numbers of nodes (8-CQUAD8, 9-CQUAD) and the CQUADR and
CQUADX are axi-symmetric versions of the CQUAD4, and CQUAD8 respectively. 
However, the ``Area()``, ``Mass()``, ``Density()``, etc. methods are calculated
in the the same way for each card (although the axi-symmetric cards return mass
per unit theta).  The last thing to note is ``rawFields`` and ``reprFields`` are
very important to how the code integrates.

``rawFields`` is used to check if a duplicated card is the same as another card
and is also used for testing.  After reading and writing, reading back in,
writing back out, reading back in, if the fields are the  same, then there's
likely no error in reading a card (fields can still be missed while reading, so
it's not perfect). ``rawFields`` returns a list of the fields (similar to the 
list-esque card object from before).

``reprFields`` is analogous to the ``__repr__()`` method, and is an abbreviated
way to write the card. For example, the ``T1, T2, T3``, and ``T4`` values
(thickness at nodes 1, 2, 3, 4) are generally 0.0 and instead are set at an
elemental level using the PSHELL card.  If these fields were printed, the CQUAD4
card would be a two line card instead of a one line card.  ``reprFields`` is
used instead of ``__repr__()`` in order to be able to write the card in large
field or small field format.  Defaults are generally not written by the
``__repr__()`` method, but are written for certain fields (e.g. the ``xyz``
fields on the GRID card).

To get the CQUAD4, with an element ID of 1000, you would type::

 elem = model.elements[1000]

or::

 elem = model.Element(1000)

to use the function.

Then to print the card, type::

 print(elem)

to see the Nastran formatted card.  The ``__repr__()`` method is defined in 
``bdf/cards/baseCard.py`` the :py:class:`pyNastran.bdf.cards.baseCard` class 
(which is used by the :py:class:`pyNastran.bdf.cards.baseCard.Element` class 
also defined in ``baseCard.py``).


:mod:`shell`: Cross-Referencing the CQUAD4 Object
--------------------------------------------------

Previously, it was mentioned that the square plate model built with quads and 
triangles had a thickness and a material. The nodes of the elements also have 
positions.  The nodes further be defined in a rectangular, cylindrical, 
or spherical coordinate system, so to get the mass of an element is actually
quite involved.  Creating a function to access the mass becomes possible without
passing the entire model object around to every function through the use of 
cross-referencing.

Cross Referencing takes a CQUAD4 card and replaces the GRID references with actual
GRID cards.  The GRID cards in turn reference two COORDx (CORD1R, CORD2R, CORD1C,
COR2DC, CORD1S, CORD2S) cards, which also may reference two CORDx cards. 
The CQUAD4 references a PSHELL or PCOMP card.  The PSHELL references a single
MAT1 card, and as mentioned before the PCOMP card may reference one or more 
MAT1/MAT8 cards.  In order to calculate something simple like the mass of the
CQUAD4 requires the formula:

.. math::
 m = A \left( t \rho + \frac{nsm}{A} \right)

for a PSHELL or:

.. math::
 m = A \left( \sum_{i=0}^{i=1}{t\rho} + \frac{nsm}{A} \right)

for a PCOMP.


By using classes and functions, it's easy to just call the ``element.MassPerArea()``
method and get the proper data to apply the formula.  Similarly, the 
``element.Area()`` method calls the ``node.Position()`` method to get the node
in the global XYZ coordinate frame and can then find the area using vectors 
in a 3D space:

.. math::
 A=\frac{1}{2} | (n_1-n_3) \times (n_2-n_4) |

(see http://en.wikipedia.org/wiki/Quadrilateral).


:mod:`crossReference`: Cross-Referencing Process
-------------------------------------------------

Cross referencing must first be done on the coordinate cards.  Then, once they're
done, the nodes are cross referenced. Once this is done, the coordinate systems
may be resolved (CORD1x cards reference GRID cards).  Then elements, properties,
materials, loads, boundary conditions, aerodynamic cards, thermal, dynamic, 
etc. cards are mapped. The order doesn't matter, but CORD1x cards and GRID cards
must be mapped first.

Cross Referencing is performed by looping over the card objects and calling the
``card.crossReference()`` method.  This will setup all cross-referencing and 
a full list of the status of various cards is listed in ``bdf_crossReferencing.txt``.

:mod:`writeMesh.py`: Writing the BDF
-----------------------------------------

The BDF is written by looping through all the objects and calling the
``card.__repr__()`` method by typing ``str(card)``.

Currently, only small field format is supported for writing.  The list from
``reprFields()`` is passed into ``fieldWriter.py`` function ``printCard(listObj)``
and it dynamically figures out how to write the card based on the data type. 
For float values, the highest precision 8-character width field will be used 
even if it uses Nastran's strange syntax of "1.2345+8" to represent 
a more standard "1.2345e+08".
