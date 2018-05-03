.. _xref-test_bdf:

========
test_bdf
========

``test_bdf`` is a very useful tool to check if pyNastran will
read/write a BDF.  Additionally, ``test_bdf`` is very useful when
debugging bad decks.

It'll crash when ``INCLUDE`` files aren't found.

It'll crash if you have invalid types.
It'll tell you the card, the field name and number as well.

It'll crash if you have missing nodes, coordinate systems, elements, 
properties, materials, loads, optimization cards, etc.

It'll crash if your TRIM card that's referenced in your Case Control is
over/underdefined (by looking at the data on your AESTAT, AESURF, AELINK,
AEPARM, SUPORT1, and SUPORT cards), but not if it's unused.

Afterwards, you can then use a program like
`Beyond Compare <http://scootersoftware.com/>`_ to see the differences.
Most of the changes you'll see are related to positioning of cards,
precision of the data fields, and values that were the default value
are left blank in order to simplify the file.

**An ENDDATA card is highly recommended** as it's an easy way to tell
if all the lines were parsed properly.


Details
=======
After installing pyNastran, open a command prompt in a folder you have
a BDF you want to test and type:

.. code-block:: conosle

  >>> test_bdf fem.bdf > test.out

The file test.out will be created...

.. code-block:: conosle

  INFO:    fname=bdf.pyc                lineNo=371    ---starting BDF.read of fem.bdf---
  INFO:    fname=bdf.pyc                lineNo=589    reject_card_name = |TEMPD|
  INFO:    fname=bdf.pyc                lineNo=589    reject_card_name = |CTRIA3|
  INFO:    fname=bdf.pyc                lineNo=384    ---finished BDF.read of fem.bdf---
  INFO:    fname=write_mesh.pyc         lineNo=68     ***writing fem02.bdf_out


  INFO:    fname=bdf.pyc                lineNo=371    ---starting BDF.read of fem.bdf_out---
  INFO:    fname=bdf.pyc                lineNo=589    reject_card_name = |TEMPD|
  INFO:    fname=bdf.pyc                lineNo=384    ---finished BDF.read of fem.bdf_out---
  INFO:    fname=write_mesh.pyc         lineNo=68     ***writing fem02.bdf_out2


  diff_keys1=[] diff_keys2=[]
     key=CHEXA   value1=52   value2=52
     key=CPENTA  value1=52   value2=52
     key=ENDDATA value1=1    value2=1
     key=GRID    value1=135  value2=135
    *key=INCLUDE value1=1    value2=0
     key=MAT4    value1=1    value2=1
     key=NLPARM  value1=1    value2=1
     key=PARAM   value1=2    value2=2
     key=PSOLID  value1=1    value2=1
     key=SPC     value1=28   value2=28
     key=TEMP    value1=18   value2=18
    -key=TEMPD   value1=1    value2=1
    *key=CTRIA3  value1=1    value2=0
  ----------------------------------------------------------------



It's broken into 3 sections:
 - **Section 1** - tells you what file was read and what cards were rejected
 - **Section 2** - uses the output of Section 1 to test again (helps to
                   verifies fields weren't shifted; not a perfect test)
 - **Section 3** - prints out the cards that were found in the BDF along
                   with the number of each.
 - A star (e.g. *key=CTRIA3) indicates a cards was not written out and
   you should be careful when using the code on the example.  The CTRIA3
   card was lost.  Most likely the CTRIA3 shares an elementID with
   another element.  Note that INCLUDE files always have stars by them.
 - A dash indicates the cards was rejected and will be echoed to the
   output BDF.  Looking in Sections 1 & 2, tells us that the TEMPD card
   was ignored in both cases.


You can also test the code without using cross-referencing

.. code-block:: conosle

  >>> test_bdf -x fem.bdf > test.out

The full calling signature of ``test_bdf`` is shown below.


Duplicated Cards
================

Sometimes this will happen

.. code-block:: conosle

  diffKeys1=[] diffKeys2=[]
    *key=CROD    value1=2   value2=4

If you look at the definition of the card, you'll see you can define 2
instances of them on a single Nastran card.  pyNastran rewrites them in
long form (unless specified).

The list of duplicate defined cards include (not a full list):
 * CORD1R
 * CORD1C
 * CORD1S
 * CROD
 * DAREA
 * SPOINT
 * PELAS
 * PDAMP
 * PMASS
 * SPOINT (writes in short form)

Modified Cards
==============
Sometimes this will happen

.. code-block:: console

  diffKeys1=[] diffKeys2=[]
    *key=FREQ1    value1=2   value2=0
    *key=FREQ     value1=0   value2=1

If you look at the definition of the cards, you'll see that if two
FREQx cards have the same ID, then they will both be used during a
frequency analysis.  The simplest form of the card will be used (a FREQ
card in this case), even if it means combining and/or changing card
types.

The list of modified cards include:
 * FREQ
 * FREQ1
 * FREQ2

When things go wrong
====================

Try:

.. code-block:: console

  >>> test_bdf -x fem.bdf > test.out

If it's still failing, you probably are have:
  * an unsupported card in your deck
  * the include files are too complicated (Nastran allows include files
    to be referenced from any previously referenced directory).
    pyNastran only references from the base input file.  Stick all your
    include files in the same folder.
  * you're missing cards


Calling Signature
-----------------

.. code-block:: console

  C:\work>test_bdf --help

  Usage:
    test_bdf [-q] [-D] [-i] [-e E] [-x] [-p] [-c] [-L] [-f] BDF_FILENAME
    test_bdf [-q] [-D] [-i] [-e E] [-x] [-p] [-c] [-L] [-d] [-f] BDF_FILENAME
    test_bdf [-q] [-D] [-i] [-e E] [-x] [-p] [-c] [-L] [-l] [-f] BDF_FILENAME
    test_bdf [-q] [-D] [-i] [-e E] [-p] [-r] [-f] BDF_FILENAME
    test_bdf [-q] [-D] [-i] [-e E] [-x] [-p] [-s] [-f] BDF_FILENAME
    test_bdf -h | --help
    test_bdf -v | --version

  Positional Arguments:
    BDF_FILENAME   path to BDF/DAT/NAS file

  Options:
    -q, --quiet    prints debug messages (default=False)
    -x, --xref     disables cross-referencing and checks of the BDF.
                    (default=False -> on)
    -p, --punch    disables reading the executive and case control decks in the BDF
                   (default=False -> reads entire deck)
    -c, --check    disables BDF checks.  Checks run the methods on
                   every element/property to test them.  May fails if a
                   card is fully not supported (default=False)
    -l, --large    writes the BDF in large field, single precision format (default=False)
    -d, --double   writes the BDF in large field, double precision format (default=False)
    -L, --loads    Disables forces/moments summation for the different subcases (default=False)
    -r, --reject   rejects all cards with the appropriate values applied (default=False)
    -D, --dumplines  Writes the BDF exactly as read with the INCLUDES processed (pyNastran_dump.bdf)
    -i, --dictsort  Writes the BDF with exactly as read with the INCLUDES processed (pyNastran_dict.bdf)
    -f, --profile   Profiles the code (default=False)
    -s, --stop      Stop after first read/write (default=False)
    -e E, --nerrors E  Allow for cross-reference errors (default=100)
    -h, --help     show this help message and exit

What test_bdf does
------------------
``test_bdf`` is useful for quickly finding out if your deck has a
cross-reference error, which may not be caught by Nastran (often
leading to segfaults).  It's also useful in testing to see if your
specific set of cards handle large/double precison in Nastran (and
that the code doesn't have a bug).

The standard test is:

.. code-block:: conosle

  >>> test_bdf model.bdf

which will:
 1. read the model

 2. card checks

    - types of each card field
    - validate values of strings (e.g. HAT is a valid PBARL)
    - check range of certain values (e.g. CONM2 must be positive semidefinite)

 3. cross-reference the model
   - Passing requires all traced cards exist:

     - Structural:

       - GRID -> Coord
       - Coord -> Coord
       - Element -> Property
       - Property -> Material
       - Material -> Nonlinear Material
       - Load -> Load cards
       - DLOAD -> DLoad cards

     - Optimization

       - DRESP -> DRESP, DESVAR
       - DEQATN -> DRESP, DESVAR

     - Aero:

       - CAERO -> PAERO
       - SPLINE -> CAERO, SET

 4. validate the methods for the various cards

   - Area
   - Mass
   - DEQATN format

 5. write the model

 6. read it back in

 7. compare it to the original

 8. calculate mass/forces/moments for all the subcases

