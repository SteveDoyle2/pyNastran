# Future Development #

If you're not familar with Nastran, reading "[What is Nastran](http://code.google.com/p/pynastran/wiki/What_Is_Nastran)" will help explain the importance of the software.  Also, if you don't know Python, don't worry, it's very easy to pick up.

Major areas of future development include:

| **Name**| **Knowledge of Python** | **Knowledge of C++** | **Level of Difficulty** |
|:|:------------------------|:---------------------|:------------------------|
| BDF Reader/Writer    | low    | None   | low    |
| F06 Reader           | low    | None   | low    |
| F06 Writer           | low    | None   | low    |
| OP2 Writer           | medium | None   | medium |
| OP2 Reader           | medium | None   | medium |
| FEA Documentation    | None   | None   | high   |
| BDF Reader in C++    | low    | medium | medium |
| GUI Development      | medium | None   | medium |
| Code Aster Converter | low    | None   | high   |
| Code Aster Parser    | high   | None   | medium |
| Code Documentation   | low    | None   | low    |
| Lines of Code Reduction  | low    | None   | low    |

# BDF Reader/Writer #
There are hundreds of cards supported by Nastran.  The cards follow a simple format, and adding new cards is easy.  The biggest focus here is testing and there are about 3000 models to use to compare to.

# Software Test Development #
Large sections of the code are fixed and are unlikely to change much going forward and don't need software testing as rigorously.  However, there are a few areas of the code that are often modified and a small error can be difficult to catch.  Additional software testing would help decrease end of version debugging.

# F06 Writer #
This is an easy task, but is more tedious than other tasks as it requires one to compare files to see if there are any differences.  It is important that the spacing of the data is correct in order to plug-in to commercial software properly.  A program like Beyond Compare is very useful.

# OP2 Reader #
Under the hood the binary OP2 Reader does a lot with Python's object oriented capability.  A dictionary containing all the important information (dataCode) is stored.  Many more result types need to be added and some major rewrites to the code are required.  Binary file reading is required, but is easy to learn thanks to the helpful tools and format documentation.

# OP2 Writer #
The OP2 Writer interfaces with the data the OP2 Reader has already read.  Presumably, a user will populate these objects with results from a code like Code Aster and then write it to the Nastran OP2.  Familiarity with the OP2 Reader Objects will be necessary in doing this task.  However, they're fairly straight-forward once you learn the few quirks they all have.    The biggest focus here is testing and there are about 3000 models to use to compare to.

# Code Aster Converter #
There are many Nastran solutions and many Code Aster solutions.  Though the cards are not necessarily one-to-one, the converter will enable running a widely used, open-source FEA code.  In the places there are not one-to-one, user judgement will be required.  Most of the documentation is in French, which would make it difficult for many users.  Access to Linux (potentially through a virtual machine using a program like [VMWare Player](http://www.vmware.com/products/player/)) is necessary.  For this task, understanding of programming is less important than being able to read the manuals, creating simple problems, running Code Aster, and trying to match the results.

# Code Aster Results Parser #
The various results from Code Aster need to be parsed and subsequently written to the OP2.  For the most part, there should be a one-to-one correlation between displacement in Code Aster and Nastran.  In the places there are not, user judgement will be required.  Access to Linux (potentially through a virtual machine using a program like [VMWare Player](http://www.vmware.com/products/player/)) is necessary.  Having programmed in Python (or any other language) before will help to allow one to develop a large scale parser with less difficulty.

# GUI Development #
The existing GUI is incapable of displaying nodal and centroidal results without corrupting the data objects.  Fixing this would greatly improve the capability of the software.  Additionally, marker (arrow) plots, deflection plots, 2D plots, element/node selection and other general GUI enhancements would greatly improve the usefulness of the program.  The output parsing capability already exists.  Currently, the code is written in wxPython and VTK.  However, switching wxPython out for PyQt4 may be a better choice.

# Documentation #
Many of the functions in the various classes are undocumented.  The HTML documentation would be more useful if every parameter was documented.  Luckily, most of the code is pretty straight-forward and the areas that aren't are more likely to be documented aren't as complicated.  Areas that a user would interface with are higher priority.  Additionally, documentation of the majority of parameter is available in available PDF files.

# Lines of Code Reduction #
The current development version of pyNastran has 52k lines of code spanning 222 files.  Reducing this by 10% would help tremendously.

Some statistics:
| **Line Count**  | **Number of Lines** |
|:----------------|:--------------------|
| Total lines   | 72150 |
| Code lines    | 52101 |
| Comment lines | 10831 |
| Blank lines   | 9218  |