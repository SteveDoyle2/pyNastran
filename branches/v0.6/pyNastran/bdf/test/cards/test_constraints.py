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
from pyNastran.bdf.bdf import BDFCard, SUPORT, SUPORT1

class TestConstraints(unittest.TestCase):
    def test_suport1_01(self):
        card = ['SUPORT1', '1']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT1(card)

        card = ['SUPORT1', '1', '2']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')

        card = ['SUPORT1', '1', '2', '432']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')

        card = ['SUPORT1', '1', '2', '432', 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)

        card = ['SUPORT1', '1', '2', None, 3]
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')

        card = ['SUPORT1', '1', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT1(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')

    def test_suport_01(self):
        card = ['SUPORT']
        card = BDFCard(card)
        with self.assertRaises(AssertionError):  # too short
            SUPORT(card)

        card = ['SUPORT', '2']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.Cs[0], '234')
        self.assertEqual(len(con.IDs), 1)
        self.assertEqual(len(con.Cs), 1)

        card = ['SUPORT', '2', '432', 3]
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', None, 3]
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[0], 2)
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[0], '0')
        self.assertEqual(con.Cs[1], '0')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

        card = ['SUPORT', '2', '432', 3, '1325']
        card = BDFCard(card)
        con = SUPORT(card) # default
        self.assertEqual(con.IDs[1], 3)
        self.assertEqual(con.Cs[1], '1235')
        self.assertEqual(len(con.IDs), 2)
        self.assertEqual(len(con.Cs), 2)

if __name__ == '__main__':
    unittest.main()