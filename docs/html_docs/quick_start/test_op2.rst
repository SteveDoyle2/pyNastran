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
