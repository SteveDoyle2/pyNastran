# Introduction #

In order to check the software will work read/write a BDF with a given set of cards, it's useful to run the test\_bdf utility.  You can then use a program like "Beyond Compare" to see the differences for yourself.  Most of the changes you'll see are related to positioning of cards, precision of the data fields, and values that were the default value are left blank in order to simplify the file.


**An ENDDATA card is highly recommended** as it's an easy way to tell if all the lines were parsed properly.


# Details #

After installing pyNastran, open a command prompt in a folder you have a BDF you want to test and type:

```

>> test_bdf cylinder.bdf > test.out

```

The file test.out will be created...

```
INFO:    fname=bdf.pyc                lineNo=371    ---starting BDF.read of cylinder.bdf---
INFO:    fname=bdf.pyc                lineNo=589    RejectCardName = |TEMPD|
INFO:    fname=bdf.pyc                lineNo=589    RejectCardName = |CTRIA3|
INFO:    fname=bdf.pyc                lineNo=384    ---finished BDF.read of cylinder.bdf---
INFO:    fname=bdf_writeMesh.pyc      lineNo=68     ***writing cylinder02.bdf_out


INFO:    fname=bdf.pyc                lineNo=371    ---starting BDF.read of cylinder.bdf_out---
INFO:    fname=bdf.pyc                lineNo=589    RejectCardName = |TEMPD|
INFO:    fname=bdf.pyc                lineNo=384    ---finished BDF.read of cylinder.bdf_out---
INFO:    fname=bdf_writeMesh.pyc      lineNo=68     ***writing cylinder02.bdf_out2


diffKeys1=[] diffKeys2=[]
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
```


It's broken into 3 sections:
  * **Section 1** - tells you what file was read and what cards were rejected.
  * **Section 2** - uses the output of Section 1 to test again (helps to verifies fields weren't shifted; not a perfect test)
  * **Section 3** - prints out the cards that were found in the BDF along with the number of each.
  * A star indicates a cards was not written out and you should be careful of using the code on the example.  The CTRIA3 card was lost.  Most likely the CTRIA3 shares an elementID with another element.  Note that INCLUDE files always have stars by them.
  * A dash indicates the cards was rejected and will be echoed to the output BDF.  Looking in Sections 1 & 2, tells us that the TEMPD card was ignored in both cases.


You can also test the code without using cross-Referencing
```
test_bdf -x fem.bdf > test.out
```

```
The full list of arguments is:

usage: test_bdf [-h] [-q] [-x] [-v] bdfFileName

Tests to see if a BDF will work with pyNastran.

positional arguments:
  bdfFileName    path to BDF/DAT file

optional arguments:
  -h, --help     show this help message and exit
  -q, --quiet    Prints debug messages (default=False)
  -x, --xref     Disables cross-referencing of the BDF
  -v, --version  show program's version number and exit
```

## Duplicated Cards ##

Sometimes this will happen
```
diffKeys1=[] diffKeys2=[]
  *key=CROD    value1=2   value2=4
```
If you look at the definition of the card, you'll see you can define 2 instances of them on a single Nastran card.  pyNastran rewrites them in long form (unless specified).

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

## Modified Cards ##
Sometimes this will happen
```
diffKeys1=[] diffKeys2=[]
  *key=FREQ1    value1=2   value2=0
  *key=FREQ     value1=0   value2=1
```
If you look at the definition of the cards, you'll see that if two FREQx cards have the same ID, then they will both be used during a frequency analysis.  The simplest form of the card will be used (a FREQ card in this case), even if it means combining and/or changing card types.

The list of modified cards include:
  * FREQ
  * FREQ1
  * FREQ2

## When things go Wrong ##
Try:
```
test_bdf -x fem.bdf > test.out
```
If it's still failing, you probably are have:
  * an unsupported card in your deck
  * the include files are too complicated (Nastran allows include files to be referenced from any previously referenced directory).  pyNastran only references from the base input file.  Stick all your include files in the same folder.
  * you're missing cards