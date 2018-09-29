==================
Installed Programs
==================

--------------------------
List of Installed Programs
--------------------------

Various Command line utilities are installed by pyNastran:

 - pyNastranGUI_
 - test_bdf_
 - test_op2_
 - test_op4_
 - bdf_
 - f06_
 - format_converter_

------------
pyNastranGUI
------------

The Graphical User Interface (GUI) looks like:

.. image:: ../../../pyNastran/gui/images/qt.png

.. code-block:: console

See :doc:`gui` for more details.

--------
test_bdf
--------
Runs through various checks on a BDF that Nastran doesn't do.  Verifies your model is referenced properly.  Creates a summary table.

See :doc:`test_bdf` for more details.

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

--------
test_op2
--------
Runs through various checks on an OP2 file.  Creates a summary table.

See :doc:`test_op2` for more details.

.. code-block:: console

  >>> test_op2 -c ISat_Dploy_Sm.op2

  ...
  
  ---stats for isat_dploy_sm.op2---
  eigenvectors[1]
    isubcase = 1
    type=RealEigenvectorArray ntimes=203 nnodes=5367, table_name=BOPHIG
    data: [t1, t2, t3, r1, r2, r3] shape=[203, 5367, 6] dtype=float32
    gridTypes
    sort1
    modes = [  1   2   3 ..., 201 202 203]
    eigrs = [        0.         0.         0. ...,  11912279.  12843625.  13110797.]
    mode_cycles = [   0.       0.       0.    ...,  549.31   570.38   576.282]

  cbar_force[1]
    type=RealCBarForceArray ntimes=203 nelements=790
    data: [ntimes, nnodes, 8] where 8=[bending_moment_a1, bending_moment_a2, bending_moment_b1, bending_moment_b2, shear1, shear2, axial, torq
  ue]
    data.shape = (203, 790, 8)
    element name: CBAR-34
    sort1
    modes = [  1   2   3 ..., 201 202 203]
    eigrs = [        0.         0.         0. ...,  11912279.  12843625.  13110797.]
    cycles = [   0.       0.       0.    ...,  549.31   570.38   576.282]

  ctria3_stress[1]
    type=RealPlateStressArray ntimes=203 nelements=32 nnodes_per_element=1 nlayers=2 ntotal=64
    data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
    data.shape=(203L, 64L, 8L)
    element type: CTRIA3
    s_code: 1
    sort1
    modes = [  1   2   3 ..., 201 202 203]
    eigrs = [        0.         0.         0. ...,  11912279.  12843625.  13110797.]
    mode2s = [0 0 0 ..., 0 0 0]
    cycles = [   0.       0.       0.    ...,  549.31   570.38   576.282]

  cquad4_stress[1]
    type=RealPlateStressArray ntimes=203 nelements=4580 nnodes_per_element=1 nlayers=2 ntotal=9160
    data: [ntimes, ntotal, 8] where 8=[fiber_distance, oxx, oyy, txy, angle, omax, omin, von_mises]
    data.shape=(203L, 9160L, 8L)
    element type: CQUAD4
    s_code: 1
    sort1
    modes = [  1   2   3 ..., 201 202 203]
    eigrs = [        0.         0.         0. ...,  11912279.  12843625.  13110797.]
    mode2s = [0 0 0 ..., 0 0 0]
    cycles = [   0.       0.       0.    ...,  549.31   570.38   576.282]

  eigenvalues[ISAT_SM_DEPLOYED MODES TO 400 HZ]
    type=RealEigenvalues neigenvalues=203
    title, extraction_order, eigenvalues, radians, cycles, generalized_mass, generalized_stiffness

Or more simply:

.. code-block:: console

  >>> test_op2 -ct ISat_Dploy_Sm.op2

  ---stats for isat_dploy_sm.op2---
  eigenvectors[1]
  cbar_force[1]
  ctria3_stress[1]
  cquad4_stress[1]
  eigenvalues[u'ISAT_SM_DEPLOYED MODES TO 400 HZ']

--------
test_op4
--------
Limited checker for testing to see if an OP4 file will load.

.. code-block:: console

 >>> test_op4 --help
 Usage:
 test_op4 [-q] [-o] OP4_FILENAME
   test_op4 -h | --help
   test_op4 -v | --version

 Tests to see if an OP4 will work with pyNastran

 Positional Arguments:
   OP4_FILENAME         Path to OP4 file

 Options:
   -q, --quiet          Suppresses debug messages (default=False)
   -o, --write_op4      Writes the op2 to fem.test_op4.op4 (default=True)
   -h, --help           Show this help message and exit
   -v, --version        Show program's version number and exit

---
bdf
---

Interface to various BDF-related command line tools

.. code-block:: console

  >>> bdf --help

  Usage:
    bdf merge         (IN_BDF_FILENAMES)... [-o OUT_BDF_FILENAME]
    bdf equivalence   IN_BDF_FILENAME EQ_TOL
    bdf renumber      IN_BDF_FILENAME [-o OUT_BDF_FILENAME]
    bdf mirror        IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [--plane PLANE] [--tol TOL]
    bdf export_mcids  IN_BDF_FILENAME [-o OUT_CSV_FILENAME] [--no_x] [--no_y]
    bdf split_cbars_by_pin_flags  IN_BDF_FILENAME [-o OUT_BDF_FILENAME] [-p PIN_FLAGS_CSV_FILENAME]
    bdf bin          IN_BDF_FILENAME AXIS1 AXIS2 [--cid CID] [--step SIZE]

    bdf merge         -h | --help
    bdf equivalence   -h | --help
    bdf renumber      -h | --help
    bdf mirror        -h | --help
    bdf export_mcids  -h | --help
    bdf split_cbars_by_pin_flags  -h | --help
    bdf bin          -h | --help
    bdf -v | --version

---
f06
---

Interface to various F06-related command line tools

.. code-block:: console

  >>> f06 --help

  Usage:
    f06 plot_145 F06_FILENAME [--noline] [--modes MODES] [--subcases SUB] [--xlim FREQ] [--ylim DAMP]

    f06 plot_145 -h | --help
    f06 -v | --version

----------------
format_converter
----------------
Converts between various common formats, typically using Nastran as a common format.
This allows methods like nodal equivalencing to be written once.

.. code-block:: console

  >>> format_converter --help

  Usage:
    format_converter nastran <INPUT> <format2> <OUTPUT> [-o <OP2>]
    format_converter <format1> <INPUT> tecplot <OUTPUT> [-r RESTYPE...] [-b] [--block] [-x <X>] [-y <Y>] [-z <Z>]
    format_converter <format1> <INPUT> stl     <OUTPUT> [-b]
    format_converter <format1> <INPUT> <format2> <OUTPUT>
    format_converter -h | --help
    format_converter -v | --version

  Options:
    format1        format type (nastran, cart3d, stl, ugrid, tecplot)
    format2        format type (nastran, cart3d, stl, ugrid, tecplot)
    INPUT          path to input file
    OUTPUT         path to output file
    -o OP2, --op2 OP2  path to results file (nastran-specific)
                   only used for Tecplot (not supported)
    -x X, --xx X   Creates a constant x slice; keeps points < X
    -y Y, --yy Y   Creates a constant y slice; keeps points < Y
    -z Z, --zz Z   Creates a constant z slice; keeps points < Z
    --block        Writes the data in BLOCK (vs. POINT) format
    -r, --results  Specifies the results to write to limit output
    -b, --binary   writes the STL in binary (not supported for Tecplot)
    -h, --help     show this help message and exit
    -v, --version  show program's version number and exit

  Notes:
    Nastran->Tecplot assumes sequential nodes and consistent types (shell/solid)
    STL/Tecplot supports globbing as the input filename
    Tecplot slicing doesn't support multiple slice values and will give bad results (not crash)
    UGRID outfiles must be of the form model.b8.ugrid, where b8, b4, lb8, lb4 are valid choices and periods are important

Example:

.. code-block:: console

    >>> format_converter tecplot tecplot.*.plt tecplot.tecplot_joined.plt -x 0.0 -y 0.0 -z 0.0
    >>> format_converter nastran fem.bdf stl fem.stl -b
    >>> format_converter nastran fem.bdf cart3d fem.tri
    >>> format_converter stl model.*.stl nastran fem.bdf
