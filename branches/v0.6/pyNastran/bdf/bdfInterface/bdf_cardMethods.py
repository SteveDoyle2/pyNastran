## GNU Lesser General Public License
## 
## Program pyNastran - a python interface to NASTRAN files
## Copyright (C) 2011-2012  Steven Doyle, Al Danial
## 
## Authors and copyright holders of pyNastran
## Steven Doyle <mesheb82@gmail.com>
## Al Danial    <al.danial@gmail.com>
## 
## This file is part of pyNastran.
## 
## pyNastran is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## pyNastran is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License
## along with pyNastran.  If not, see <http://www.gnu.org/licenses/>.
## 
# pylint: disable=C0103,R0911,R0912,R0914,R0915
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import warnings


def isLargeField(card):
    """
    :returns: True if the card is in 16-character width fields
    """
    star_fields = ['*' in field for field in card]
    if any(star_fields):
        return True
    return False


def parse_csv(sline):
    """
    Parses a Nastran single-precision formatted CSV line
    """
    slineA = sline.split(',')[:9]
    sline2 = [''] * 9

    for (i, s) in enumerate(slineA):
        sline2[i] = s
    #sline = sline2
    #print "sline2 = ",sline2
    #sline = sline2
    #sline = sline.split(',')[0:9]  # doesnt fill all fields on line
    return sline2


def make_single_streamed_card(log, card, debug=False):
    """
    Takes a card that has been split b/c it's a multiline card
    and gets rid of the required blanks in it.
    """
    cardOut = []
    n = 0
    if debug:
        log.debug("card = %s" % card)
    for field in card:
        if n - 9 == 0:
            pass
        elif n - 10 == 0:
            n -= 10
        else:
            cardOut.append(field)
        n += 1
    #print "cardOut = ",cardOut
    return cardOut


def nastran_split(log, line, is_large_field, debug=False):
    """
    Splits a single BDF line into large field or small field format

    :param line:           the BDF line
    :type line:            str
    :param is_large_field: a flag indicating small/large field format
                           (True/False)
    :type is_large_field:  bool
    :param debug: extra developer debug
    :type debug:  bool

    :returns fields: the 9 (small) or 5 (large) fields for the line

    .. note:: CSV Format is handled by parse_csv
    .. note:: tabs are handled prior to running this
    """
    if debug:
        print("is_large_field = %s" % (is_large_field))
    if is_large_field:
        fields = [line[0:8], line[8:24], line[24:40], line[40:56], line[56:72]]
    else:  # small field
        fields = [line[0:8], line[8:16], line[16:24], line[24:32],
                  line[32:40], line[40:48], line[48:56], line[56:64],
                  line[64:72]]
    #if debug:
    #    print("  fields = ", collapse(fields))

    fields2 = []
    for (i, rawField) in enumerate(fields):
        field = rawField.strip()
        #if debug:
        #if (i==9) and field=='' or field=='*' or field=='+':
            #print "skipping * or + or empty field"
        #    pass
        #else:
        if debug:
            log.debug("i=%s rawField=|%s| field=|%s|" % (i, rawField, field))
        fields2.append(field)
    return fields2