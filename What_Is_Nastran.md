# What is Nastran #

## Brief History ##
[Nastran](http://en.wikipedia.org/wiki/Nastran) (NASA Structrual ANalysis) is a series of commercial software products originally developed under a [NASA](http://en.wikipedia.org/wiki/NASA) contract in the late 1960s by [MSC Software Corporation](http://www.mscsoftware.com/) using the [Fortran](http://en.wikipedia.org/wiki/Fortran) programming language.  It uses the [Finite Element Method](http://en.wikipedia.org/wiki/Finite_element_method) which discritizes geometry into small elements and solves large [sparse matrices](http://en.wikipedia.org/wiki/Sparse_matrix) using linear algebra to find quantities like displacement and stress in order to design [structures](http://en.wikipedia.org/wiki/Structural_engineering). It became the industry standard program in part due to MSC buying competitors and incorporating their advances into their products.  After an antitrust settlement with the [FTC](http://en.wikipedia.org/wiki/Federal_Trade_Commission) in 2002, their source code was released to various organizations.  Alternative Nastran versions were soon created including [NEi-Nastran](http://www.nenastran.com/nei-nastran) and [NX Nastran](http://www.plm.automation.siemens.com/en_us/products/velocity/femap/nxNastran/).

Additionally, other competitors have emerged including [Abaqus](http://www.3ds.com/products/simulia/portfolio/abaqus/latest-release/),  [Ansys](http://www.ansys.com/Industries/Aerospace+&+Defense/Aircraft),  [PERMAS](http://www.intes.de), and [Code\_Aster](http://www.code-aster.org/V2/spip.php?rubrique18).  Most Nastran competitors also have translators that support the BDF/OP2 file format.  PERMAS even supports the Nastran format directly.


## A Finite Element Model (FEM) ##
![http://www.smartcae.com/uploads/u4/upload/breading_1.gif](http://www.smartcae.com/uploads/u4/upload/breading_1.gif)

As the model is loaded, the displacement increases.  The contour plot shows [stress](http://en.wikipedia.org/wiki/Stress_(mechanics)) (Force/Area) on the part.  At high stresses, the part will fail.  The stress vs [strain](http://en.wikipedia.org/wiki/Deformation_(mechanics)) (change in length/initial length) [curve](http://en.wikipedia.org/wiki/Stress-strain_curve) has a linear portion and a nonlinear portion.  Structures are designed to be operated linearly (to extend the life of the part), but under extreme cases take a additional load such that there is no catastrophic failure.

## The Problem with Nastran ##
Due to the Fortran legacy, the file formats of this program are often criticized as being cubersome.  For aerospace, mechanical, and civil engineers who use this program on a daily basis to analyze aircraft, spacecraft, buildings, and cars, this is a huge burden.

## File Types ##
| **File Extension** | **Acronym** | **Contents** | **Input/Output** | **Format** |
|:-------------------|:------------|:-------------|:-----------------|:-----------|
| BDF/DAT/NAS | Bulk Data File / Data / Nastran | Geometry/Loads/Constraints | Input  | ASCII  |
| F04 | File #4        | Solution Information  | Output | ASCII  |
| F06 | File #6        | Displacement/Stress/etc.   | Output | ASCII  |
| PCH | Punch        | See Note   | Output | ASCII  |
| OP2 | Output #2      | Displacement/Stress/etc.   | Output | Binary |
| OP4 | Output #4      | Matrices and Vectors       | Input/Output | ASCII/Binary |

**Note** the PCH file may be the output of an optimization run (and uses valid BDF syntax) and also may be another version of the OP2/F06 file.  This output looks very different than any of the other formats.

Nastran files aren't very descriptive with their names...

## Explanation of Files ##
The BDF (Bulk Data File; [simple example](http://code.google.com/p/pynastran/source/browse/trunk/models/solidBending.bdf) [much harder example](http://code.google.com/p/pynastran/source/browse/trunk/pyNastran/bdf/test/unit/testA.bdf)) is the input file and defines nodes, which are connected by elements.  Loads and constraints are then applied to models.  Displacements and element stresses may then be solved for.  The BDF contains roughly 700 "cards", which are often 8-character fixed width and values like "1.0", but may also be in scientific notation "1.0E+07" or to save space "1+7".  They may also be written in 16-character fixed with or comma separated value (CSV) format.  These cards are highly interlinked and can be used to find the weight of a structure.  They are also very prone to error, but allow very complicated problems to be solved, often with the same cards.

Problems are solved using a variety of linear algebra techniques using a variety of solutions (e.g. statics, dynamics, modal, thermal, aeroelasticity, acoustic, optimization, etc.).  These matrices may be accessed from the [OP4](http://code.google.com/p/pynastran/source/browse/trunk/pyNastran/op4/test/mat_t_s1.op4) which may be an ASCII or binary file.  Using the DMAP Programming language allows a user to "jump" into the code and alter the behavior of it.  They can then define their own solution type.  The F04 file is very useful when using DMAP or optimization to debug a result.

Finally, there are two major output formats, the [F06](http://femci.gsfc.nasa.gov/sine_vib/sine_test.f06) and the OP2.  The F06 is an ASCII file and adequate for small problems and many potential users of pyNastran have written parsers for it.  However, it has a few problems that the OP2 does not.  It leaves useful data that aide in parsing output such as what solution type was run, and due to an 80 character width, many float values are significantly truncated.  The OP2 file is binary and does not have these problems.  The file size also roughly 5 times smaller for an equivalent output and OP2 parsers are much faster than F06 parsers.  However, the OP2 introduces it's own challenges by being written in a Fortran formatted block and embeds buffer tables periodically into the file.

Similar to the BDF, there are new output result tables for every solution type.  Results such as displacements, velocity, acceleration, stresses, temperatures, heat flux, element energy, element forces, grid point forces, eigenvalues, eigenvectors, failure indices (for composites).  Additionally, these tables has a corresponding complex table and random version of the table.  Reading the output of Nastran is a huge task.

# What's Really so Useful about pyNastran? #
For input file reading, pyNastran provides BDF interlinking and access to the user to modify/remove/add data to try and solve the problem they're interested in.

In terms of results processing, the number of tables in Nastran is roughly nResultTypes\*real\_complex\_random\*nResultTables = 12\*3\*200=7200.  This is an excessive number of tables, but by combining various table types into a single result (e.g. a static, a time domain and frequency domain solution) through use of classes and dictionaries, the required number of tables to parse is drastically reduced.  Finally, many tables are similar and identifying these tables helps to reduce the workload even further.