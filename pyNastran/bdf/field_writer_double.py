"""
Defines functions for double precision 16 character field writing.
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from six.moves import range

import sys
from pyNastran.utils import integer_types
from pyNastran.bdf.cards.utils import wipe_empty_fields_typed

def print_scientific_double(value):
    """
    Prints a value in 16-character scientific double precision.

    Scientific Notation:                   5.0E+1
    Double Precision Scientific Notation:  5.0D+1
    """
    if value < 0:
        Format = "%16.9e"
    else:
        Format = "%16.10e"

    svalue = Format % value
    field = svalue.replace('e', 'D')

    if field == '-0.0000000000D+00':
        field = '0.0000000000D+00'
    #assert len(field) == 16, ('value=%r field=%r is not 16 characters '
    #                          'long, its %s' % (value, field, len(field)))
    return field


def print_field_double(value):
    """
    Prints a 16-character width field

    :param value:   the value to print
    :returns field: an 16-character string
    """
    if isinstance(value, integer_types):
        field = "%16s" % value
    elif isinstance(value, float):
        field = print_scientific_double(value)
    elif value is None:
        field = "                "
    else:
        field = "%16s" % value
    if len(field) != 16:
        msg = 'field=%r is not 16 characters long...rawValue=%r' % (field, value)
        raise RuntimeError(msg)
    return field


def print_card_double(fields, wipe_fields=True):
    """
    Prints a nastran-style card with 16-character width fields.

    :param fields: all the fields in the BDF card (no trailing Nones)
    :param wipe_fields:  some cards (e.g. PBEAM) have ending fields
                         that need to be there, others cannot have them.
    .. note:: An internal field value of None or '' will be treated as
              a blank field
    .. note:: A large field format follows the  8-16-16-16-16-8 = 80
              format where the first 8 is the card name or
              blank (continuation).  The last 8-character field indicates
              an optional continuation, but because it's a left-justified
              unneccessary field, print_card doesnt use it.

    .. code-block:: python

      >>> fields = ['DUMMY', 1, 2, 3, None, 4, 5, 6, 7, 8.]
      >>> print_card_double(fields)
      DUMMY*                 1               2               3
      *                      4               5               6               7
      *       8.0000000000D+00
      *
    """
    if wipe_fields:
        fields = wipe_empty_fields_typed(fields)
    nfields_main = len(fields) - 1  # chop off the card name
    nbdf_lines = nfields_main // 8
    if nfields_main % 8 != 0:
        nbdf_lines += 1
        nextra_fields = 8 * nbdf_lines - nfields_main
        fields += [None] * nextra_fields

    try:
        out = '%-8s' % (fields[0] + '*')
    except:
        print("ERROR!  fields=%s" % fields)
        sys.stdout.flush()
        raise

    for i in range(1, len(fields)):
        field = fields[i]
        try:
            out += print_field_double(field)
        except:
            print("bad fields = %s" % fields)
            raise
        if i % 4 == 0:  # allow 1+4 fields per line
            out = out.rstrip(' ')
            if out[-1] == '\n':  # empty line
                out += '*'
            out += '\n*       '
    out = out.rstrip(' *')  # removes one continuation star
    if not out.endswith('\n'):
        out += '\n'
    return out
