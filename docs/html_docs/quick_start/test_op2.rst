.. _xref-test_op2:

========
test_op2
========

``test_op2`` verifies that the OP2 is read properly.  It's mainly a
developer debugging script as it runs the OP2 twice (with different
routines to make sure the answers are the same), but ``test_op2``
is very useful for understanding what is inside an OP2.

In general, it's recommended that in general you call ``test_op2`` like:

.. code-block:: console

  >>> test_op2 -c fem.op2

If you want an F06 file:

.. code-block:: console

  >>> test_op2 -cf fem.op2

You can skip results to minimize memory usage.  Skipping stress and rod_strain:

.. code-block:: console

  >>> test_op2 -cf fem.op2 -x stress -x rod_strain

You may also skip specific subcases (read subcases 1, 5):

.. code-block:: console

  >>> test_op2 -cf fem.op2 -s 1_5


Finally, you can extract the geometry and write a BDF.

.. code-block:: console

  >>> test_op2 -c -gn fem.op2

or

  >>> test_op2 -c --geometry --write_bdf fem.op2


Example
-------

Here, we'll determine what all the tables in the OP2 are as well as a summary of the objects

.. code-block:: console

  >>> test_op2 -c ISat_Dploy_Sm.op2

  OP2_FILENAME = 'isat_dploy_sm.op2'
  DEBUG:     fname=op2_scalar.py             lineNo=521    set_subcases - subcases = []
  DEBUG:     fname=op2_scalar.py             lineNo=521    set_subcases - subcases = []
  DEBUG:     fname=op2.py                    lineNo=385    combine=True
  DEBUG:     fname=op2.py                    lineNo=386    -------- reading op2 with read_mode=1 (array sizing) --------
  INFO:      fname=op2_scalar.py             lineNo=1398   op2_filename = 'isat_dploy_sm.op2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='PVT0'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CSTM'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=2241   table_name = 'GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='EPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='MPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM4'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BGPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='DYNAMICS'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CASECC'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='LAMA'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BOPHIG'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OES1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OEF1'
  DEBUG:     fname=op2.py                    lineNo=398    -------- reading op2 with read_mode=2 (array filling) --------
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='PVT0'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CSTM'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=2241   table_name = 'GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='EPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='MPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM4'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BGPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='DYNAMICS'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CASECC'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='LAMA'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BOPHIG'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OES1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OEF1'
  DEBUG:     fname=op2.py                    lineNo=534    combine_results
  DEBUG:     fname=op2.py                    lineNo=406    finished reading op2
  DEBUG:     fname=op2.py                    lineNo=385    combine=True
  DEBUG:     fname=op2.py                    lineNo=386    -------- reading op2 with read_mode=1 (array sizing) --------
  INFO:      fname=op2_scalar.py             lineNo=1398   op2_filename = 'isat_dploy_sm.op2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='PVT0'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CSTM'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=2241   table_name = 'GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='EPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='MPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM4'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BGPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='DYNAMICS'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CASECC'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='LAMA'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BOPHIG'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OES1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OEF1'
  DEBUG:     fname=op2.py                    lineNo=398    -------- reading op2 with read_mode=2 (array filling) --------
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='PVT0'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CSTM'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=2241   table_name = 'GPL'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='EPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='MPT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM2'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM4'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='GEOM1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BGPDT'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='DYNAMICS'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='CASECC'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='LAMA'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='BOPHIG'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OES1'
  DEBUG:     fname=op2_scalar.py             lineNo=1589     table_name='OEF1'
  DEBUG:     fname=op2.py                    lineNo=534    combine_results
  DEBUG:     fname=op2.py                    lineNo=406    finished reading op2


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


  INFO:      fname=op2.py                    lineNo=639    ---self.subcase_key---
  INFO:      fname=op2.py                    lineNo=642    subcase_id=1 : keys=[1]

Alternatively, we can call it and get a shorter summary:

.. code-block:: console

  >>> test_op2 -ct ISat_Dploy_Sm.op2

  ---stats for isat_dploy_sm.op2---
  eigenvectors[1]
  cbar_force[1]
  ctria3_stress[1]
  cquad4_stress[1]
  eigenvalues[u'ISAT_SM_DEPLOYED MODES TO 400 HZ']


Calling Signature
-----------------

.. code-block:: console

  test_op2 [-q] [-b] [-c] [-g] [-n] [-m] [-f] [-o] [-p] [-z] [-w] [-t] [-s <sub>] [-x <arg>]... OP2_FILENAME
    test_op2 -h | --help
    test_op2 -v | --version

  Tests to see if an OP2 will work with pyNastran

  Positional Arguments:
    OP2_FILENAME         Path to OP2 file

  Options:
    -b, --binarydebug     Dumps the OP2 as a readable text file
    -c, --disablecompare  Doesn't do a validation of the vectorized result
    -q, --quiet           Suppresses debug messages [default: False]
    -t, --short_stats     Short get_op2_stats printout
    -g, --geometry        Reads the OP2 for geometry, which can be written out
    -n, --write_bdf       Writes the bdf to fem.test_op2.bdf (default=False)
    -f, --write_f06       Writes the f06 to fem.test_op2.f06
    -m, --write_xlsx      Writes an XLSX to fem.test_op2.xlsx
    -o, --write_op2       Writes the op2 to fem.test_op2.op2
    -p, --profile         Profiles the code (default=False)
    -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of
                          Real/Imaginary (still stores Real/Imag); [default: False]
    -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)
    -w, --is_sort2        Sets the F06 transient to SORT2
    -x <arg>, --exclude   Exclude specific results
    -h, --help            Show this help message and exit
    -v, --version         Show program's version number and exit
