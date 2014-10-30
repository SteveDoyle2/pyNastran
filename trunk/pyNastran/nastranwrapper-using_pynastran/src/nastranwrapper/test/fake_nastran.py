#!/bin/python

import os
import sys
import shutil
import difflib

# This is to run tests without owning a copy of Nastran
# The idea is that you will run this program instead of running
# Nastran. This program will make sure that you are supplying
# the right input file and can parse the output file correctly

# The first argument will be nastran's first argument: the real input file
# The second argument will be what the input file should be
# The third argument will be what the output should be
# All arguments following are intended for the real nastran, but
# will be ignored.

def main():
    if len(sys.argv) < 5:
        print("This program takes at least 4 arguments:\n" + \
              "./fake_nastran nastran's_input ideal_input ideal_output ideal_op2")
        sys.exit(1)

    nastran_input_filename = sys.argv[1]
    ideal_input_filename = sys.argv[2]
    ideal_output_filename = sys.argv[3]
    ideal_op2_filename = sys.argv[4]

    fh = open(ideal_input_filename, "r")
    ideal_input = fh.readlines()
    fh.close()

    fh = open(nastran_input_filename, "r")
    nastran_input = fh.readlines()
    fh.close()

    # First we make sure the ideal input is the same as the nastran input
    diff = list(difflib.context_diff(ideal_input, nastran_input))
    if len(diff) > 0: # we have a problem
        print("There was a difference between the ideal input and nastran's input.")
        print()
        for line in diff:
            print(line)
        sys.exit(1)

    # Now we have to output the output files given to us in
    # the same place nastran would put it
    # First, we have to identify the directory of the input file,
    # then name the output file the same thing as the input file
    # but take off the 'bdf' and replace it with 'out' and 'op2'

    folder, filename = os.path.split(nastran_input_filename)

    # first, the 'out' aka 'f06' file
    new_output_filename = filename[:-3] + "out"
    shutil.copy(ideal_output_filename, os.path.join(folder, new_output_filename))

    # second, the 'op2' file
    new_output_filename = filename[:-3] + "op2"
    shutil.copy(ideal_op2_filename, os.path.join(folder, new_output_filename))

    # And that's that.


if __name__ == "__main__":
    main()
