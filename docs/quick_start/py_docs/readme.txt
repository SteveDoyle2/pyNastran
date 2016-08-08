Installing
----------
>>> python setup_user.py develop

List of Installed Scripts
-------------------------
pyNastranGUI
test_bdf
test_op2
test_op4
format_converter
run_nastran_double_precision


pyNastranGUI
------------
The GUI can be seen at:
  https://github.com/SteveDoyle2/pyNastran/wiki/GUI

>>> pyNastranGUI --help
Usage:
  pyNastranGUI [-f FORMAT] [-i INPUT] [-o OUTPUT...]
                  [-s SHOT] [-m MAGNIFY]
                  [-g GSCRIPT] [-p PSCRIPT]
                  [-u POINTS_FNAME...]
                  [-c][-q]
  pyNastranGUI -h | --help
  pyNastranGUI -v | --version

Options:
  -h, --help                  show this help message and exit
  -f FORMAT, --format FORMAT  format type (cart3d, lawgs, nastran, panair,
                                           plot3d, stl, tetgen, usm3d)
  -i INPUT, --input INPUT     path to input file
  -o OUTPUT, --output OUTPUT  path to output file
  -c, --console               disable HTML console output
  -g GSCRIPT, --geomscript GSCRIPT  path to geometry script file (runs before load geometry)
  -p PSCRIPT, --postscript PSCRIPT  path to post script file (runs after load geometry)
  -s SHOT, --shots SHOT       path to screenshot (only 1 for now)
  -m MAGNIFY, --magnify       how much should the resolution on a picture be magnified [default: 5]
  -u POINTS_FNAME, --user_points POINTS_FNAME               add user specified points to an alternate grid (repeatable)
  -q, --quiet                 prints debug messages (default=True)
  -v, --version               show program's version number and exit

test_bdf
--------
Runs through various checks on a BDF that Nastran doesn't do.  Verifies your model is referenced properly.  Creates a summary table.

>>> test_bdf --help
Usage:
  test_bdf [-q] [-D] [-i] [-x] [-p] [-c] [-L] [-f] BDF_FILENAME
  test_bdf [-q] [-D] [-i] [-x] [-p] [-c] [-L] [-d] [-f] BDF_FILENAME
  test_bdf [-q] [-D] [-i] [-x] [-p] [-c] [-L] [-l] [-f] BDF_FILENAME
  test_bdf [-q] [-D] [-i] [-p] [-r] [-f] BDF_FILENAME
  test_bdf [-q] [-D] [-i] [-x] [-p] [-s] [-f] BDF_FILENAME
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
  -f, --profile  Profiles the code (default=False)
  -s, --stop     Stop after first read/write (default=False)
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit


test_op2
--------
Runs through various checks on an OP2 file.  Creates a summary table.

>>> test_op2 --help
Usage:
test_op2 [-q] [-b] [-c] [-g] [-n] [-m] [-f] [-o] [-p] [-z] [-w] [-s <sub>] [-x <arg>]... OP2_FILENAME
  test_op2 -h | --help
  test_op2 -v | --version

Tests to see if an OP2 will work with pyNastran 0.8.0+dev.a941748.

Positional Arguments:
  OP2_FILENAME         Path to OP2 file

Options:
  -b, --binarydebug     Dumps the OP2 as a readable text file
  -c, --disablecompare  Doesn't do a validation of the vectorized result
  -q, --quiet           Suppresses debug messages [default: False]
  -g, --geometry        Reads the OP2 for geometry, which can be written out
  -n, --write_bdf      Writes the bdf to fem.test_op2.bdf (default=False)
  -f, --write_f06       Writes the f06 to fem.test_op2.f06
  -m, --write_xlsx      Writes an XLSX to fem.test_op2.xlsx
  -o, --write_op2       Writes the op2 to fem.test_op2.op2
  -p, --profile     Profiles the code (default=False)
  -z, --is_mag_phase    F06 Writer writes Magnitude/Phase instead of
                        Real/Imaginary (still stores Real/Imag); [default: False]
  -s <sub>, --subcase   Specify one or more subcases to parse; (e.g. 2_5)
  -w, --is_sort2        Sets the F06 transient to SORT2
  -x <arg>, --exclude   Exclude specific results
  -h, --help            Show this help message and exit
  -v, --version         Show program's version number and exit
  
test_op4
--------
Limited checker for testing to see if an OP4 file will load.

>>> test_op4 --help
Usage:
test_op4 [-q] [-o] OP4_FILENAME
  test_op4 -h | --help
  test_op4 -v | --version

Tests to see if an OP4 will work with pyNastran 0.8.0+dev.a941748.

Positional Arguments:
  OP4_FILENAME         Path to OP4 file

Options:
  -q, --quiet          Suppresses debug messages (default=False)
  -o, --write_op4      Writes the op2 to fem.test_op4.op4 (default=True)
  -h, --help           Show this help message and exit
  -v, --version        Show program's version number and exit

format_converter
----------------
Converts between various common formats, typically using Nastran as a common format.  This allows methods like nodal equivalencing to be written once.

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

Examples:
  format_converter tecplot tecplot.*.plt tecplot.tecplot_joined.plt -x 0.0 -y 0.0 z 0.0
run_nastran_double_precision
----------------------------
>>> run_nastran_double_precision --help
Bugged...no prompt

Still in development, runs:
  run_nastran_double_precision fem.bdf

with various output formats to verify there are no bugs in the write_bdf method.  Verifies it by reading the OP2.

