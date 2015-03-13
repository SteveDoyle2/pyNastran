# News

## pyNastran is Moving (3/12/2015)
Google Code is  [closing down](http://google-opensource.blogspot.com/2015/03/farewell-to-google-code.html) on January 25, 2016 and as such pyNastran is moving to [github](https://github.com/SteveDoyle2/pynastran).  New commits will be made there.

Additionally, a long awaited new release (v0.7) is coming soon.  The target is sometime withing the next month.  There are many improvements such as:
 * OP2 superelement support
 * vectorized OP2 support (uses much less memory; not implemented for Element Forces)
 * additional OP2 tables
 * improved F06 reader (the OP2 reader is still better)
 * 20 new BDF cards
 * much improved GUI with transient support (real only), a results sidebar, logging, and scripting support

As always there are many bug fixes and many new tests.

## Download Page
Google Code no longer supports new downloads.  As such, the download page will now be located at  [Sourceforge](https://sourceforge.net/projects/pynastran/files/?source=navbar).  The much improved dev version of the [GUI](http://code.google.com/p/pynastran/wiki/GUI) is located there.  v0.6.1 and earlier releases will still be located here.

## Version 0.6.1 has been released
**Version 0.6** improves BDF reading.  The reader is more robust and also requires proper BDF field formatting (e.g. a integer field can't be a float).  Additionally, cards also have a comment() method.

Marcin Gąsiorek participated in the latest pyNastran under the European Space Agency's (ESA) "Summer of Code In Space" [SOCIS](http://sophia.estec.esa.int/socis2012/?q=node/13) program.  The program provides a stipend to students to work on open-source projects.  He did a great job of simplifying code and creating nicer documentation.

Also, if anyone knows how to setup a project with readthedocs.org, we can get pyNastran documentation setup online.

## Additional Info
Note the wiki refers to the most current version.
Make sure your accessing this site through https://github.com/SteveDoyle2/pynastran.

## Questions
If you have a question about how to do something with pyNastran, create a new issue or post in the [forum](http://groups.google.com/group/pynastran-discuss) or the [developer forum](https://groups.google.com/group/pynastran-dev) for questions/discussion releated to the dev version.  Also, if you have any small test problems that have issues, please make a new issue and attach the file.

# Overview
The pyNastran software interfaces to Nastran's complicated input and output files and provides a simplified interface to read/edit/write the various files.  The software is compatible currently being used on Windows, Linux, and Mac.

The **BDF reader/editor/writer** supports about 210 cards including coordinate systems.  Card objects have methods to access data such as Mass, Area, etc.  The BDF writer writes a small field formatted file, but makes full use of the 8-character Nastran field.  The OpenMDAO BDF parametrization syntax is also supported.

The **OP2 reader** supports static/transient/frequency results, which unless analyzing random responses this won't be a problem.  It also supports *F06 Writing* for most of the objects.  Results include: displacement, velocity, acceleration, temperature, eigenvectors, eigenvalues, SPC forces, MPC forces, grid point forces, load vectors, applied loads, strain energy, as well as stress and strain.

The **F06 reader** works for simple problems, but it's still preliminary.  At this point, you should just use the OP2 reader.  It's faster, more robust, and supports more results.  The F06 reader is more used as a verification tool for the OP2 reader.

The **Cython OP4 reader** supports ASCII and binary dense matrices.  The **Python OP4 reader** supports ASCII and binary sparse and dense matrices.

A simple [GUI](http://code.google.com/p/pynastran/wiki/GUI) has been developed that can view BDF models and display static stress results from the OP2.


# Project Features
 * [BDF](http://code.google.com/p/pynastran/wiki/BDF) reading/editing/writing
 * [F06](http://code.google.com/p/pynastran/wiki/F06) reading/writing
 * [OP2](http://code.google.com/p/pynastran/wiki/OP2) reading, F06 writing
 * [OP4](http://code.google.com/p/pynastran/wiki/OP4) dense matrix reading
 * [GUI](http://code.google.com/p/pynastran/wiki/GUI) for displaying a BDF and OP2 stresses
 * About 180 supported BDF cards; some with methods (e.g. mass, area, volume, centroid, thickness, length, nsm)
 * [install](http://code.google.com/p/pynastran/wiki/InstallationGuide) script
 * automated [documentation](http://code.google.com/p/pynastran/wiki/DeveloperDocs user]  using Doxygen
 * list of <A href="https://pynastran.googlecode.com/svn/tags/v0.5.0_py27/bdf_readWrite.txt">BDF cards</A>
 * list of <A href="https://pynastran.googlecode.com/svn/tags/v0.5.0_py27/op2_reading.txt">OP2 Tables</A> 
 * list of <A href="https://pynastran.googlecode.com/svn/tags/v0.5.0_py27/f06_readWrite.txt">F06 Tables</A> 
 * more technical version of features see the <A href="https://pynastran.googlecode.com/svn/tags/v0.5.0_py27/releaseNotes.txt"> Release Notes</A> 

## GUI Features
 * Tested on Windows & Linux
 * Buttons for loading a BDF and an OP2
 * Limited results viewing (only stress)
 * Show Element Edges (using command line argument)

## BDF Features
 * 210 BDF Cards stored in "card objects"
 * "Card Objects" store relevant card data (entries are typically the same as the QRG) and when printed, automatically print in Nastran small field format (8-character fields)
 * reading/editing methods for case control decks
 * add/delete methods for subcases, nodes, elements, properties, loads, SPCs/MPCs, spoints, case control deck and more...
 * automatic checks so you don't overwrite nodes and stuff like that
 * INCLUDE file support
 * small/large/mixed field, tabs, and CSV formatted lines supported

## BDF Limitations
 * some card formats not supported (mixed tab & CSV, replication)
 * no option to write large field formatted cards
 * there's probably about 500 cards in Nastran, but the majority of the standard cards are supported

## OP2 Features
 * Multiple subcases
 * F06 Writer support for most results
 * extract results at specified times using setTransientTimes (good for large models)
 * Allows for only selected subcase results to be extracted (e.g. subcase 1 & 100, but not 2)

 * Static & Transient solution support (SOL 101, 144, 200, etc.)
   * Any result that has a SOL 101 table.  For example, displacement for buckling works, but modal displacements are not supported.  In general, if you use RESULT-TYPE(REAL/PHASE,SORT1)=ALL, you'll be fine.  SORT2 works most of the time as the OP2 is generally written in SORT1 format, but random results have problems.
 * Results (Static/Transient, Real, SORT1)
   * Displacement (real/complex)
   * Eigenvalue Order (real)
   * Eigenvectors (real/complex)
   * SPC/MPC Forces (real/complex)
   * Element Forces (real/complex)
   * Velocity (real/complex)
   * Acceleration (real/complex)
   * Temperature
   * Stress/Strain
     * Von Mises/Max Shear
     * Fiber/Curvature
   * Load Vector (real/complex)
   * Grid Point Force Balance
   * Some nonlinear stress support (CQUAD4, CTRIA3, CTETRA, CROD)

 * Linear Isotropic Elements (Stress/Strain)
   * CBAR, CBEAM, CROD, CONROD, CTUBE, CQUAD4, CTRIA3, CHEXA, CPENTA, CTETRA, CELASx, CSHEAR
 * Linear Composite Elements (Stress/Strain)
   * CQUAD4, CTRIA3
 * Nonlinear Stress
   * CQUAD4, CTRIA3, CROD, CONROD, CTUBE
 * Hyperelastic Stress
   * CQUAD4

## OP2 Limitations
 * only supports PARAM,POST,-1
 * Limited nonlinear and hyperelastic element support



If anyone makes any specific requests I'll try to incorporate them.  They need to be specific, such as read these cards from the BDF, read these results from this OP2, or write these results to an OP2.  <b>Any sample problems that you have (to test the software with) would be appreciated.  I need small examples that are comprehensive that I can add as demo problems.</b>

Check out the following if you're interested to see what's being worked on:
 * <A href="http://code.google.com/p/pynastran/source/browse/trunk/releaseNotes.txt">Release Notes</A>
