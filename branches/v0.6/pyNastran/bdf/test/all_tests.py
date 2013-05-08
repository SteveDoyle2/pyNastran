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
import unittest

from pyNastran.bdf.test.test_field_writer import TestFieldWriter
from pyNastran.bdf.test.bdf_unit_tests import TestBDF, BaseCard_Test, CaseControlTest
from pyNastran.bdf.test.cards.test_coords import *
from pyNastran.bdf.test.cards.test_constraints import *
from pyNastran.bdf.test.cards.test_rigid import *
from pyNastran.bdf.test.cards.test_beams import *
from pyNastran.bdf.test.cards.test_shells import *

if __name__ == "__main__":
    unittest.main()
