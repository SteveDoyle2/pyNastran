import unittest
from io import StringIO
import numpy as np

from cpylog import get_logger
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard, read_bdf
from pyNastran.dev.bdf_vectorized3.cards.test.utils import save_load_deck
from pyNastran.dev.bdf_vectorized3.cards.test.test_vector_shells import (
    make_dvprel_optimization,
    #make_dvcrel_optimization,
    #make_dvmrel_optimization,
)


class TestPlotElements(unittest.TestCase):
    def test_plotel_1(self):
        """tests PLOTEL"""
        model = BDF(debug=False)
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [2., 0., 0.])
        model.add_plotel(42, [5, 1], comment='plot1')
        model.add_plotel(1, [1, 10])
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def _test_plotel_2(self):
        """tests PLOTEL"""
        model = BDF(debug=False)
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [2., 0., 0.])
        model.add_card(['PLOTEL',
                        42, 5, 1,
                        45, 1, 10], 'PLOTEL', comment='plot2', ifile=None, is_list=True, has_none=True)
        #model.add_plotel(42, [5, 1], comment='plot')
        #model.add_plotel(45, [1, 10], comment='plot')
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def test_plotel3_1(self):
        """tests PLOTEL3"""
        model = BDF(debug=False)
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [1., 1., 0.])

        card_lines = ['PLOTEL3', 42,
                      5, 1, 10]
        model.add_card(card_lines, 'PLOTEL3', comment='plot1')
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def test_plotel4_1(self):
        """tests PLOTEL4"""
        model = BDF(debug=False)
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [1., 1., 0.])
        model.add_grid(11, [0., 1., 0.])

        card_lines = ['PLOTEL4', 42,
                      5, 1, 10, 11]
        model.add_card(card_lines, 'PLOTEL4', comment='plot1')
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def test_plotel6_1(self):
        """tests PLOTEL6"""
        #       10
        #     /  |
        #    /   |
        #   /    |
        # 5------1
        model = BDF(debug=False)
        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [1., 1., 0.])

        model.add_grid(15, [0.5, 0., 0.])
        model.add_grid(11, [1., 0.5, 0.])
        model.add_grid(110, [0.5, 0.5, 0.])

        card_lines = ['PLOTEL6', 42,
                      5, 1, 10,
                      15, 11, 110]
        model.add_card(card_lines, 'PLOTEL6', comment='plot1')
        model.setup(run_geom_check=True)
        save_load_deck(model)

    def test_plotel8_1(self):
        """tests PLOTEL8"""
        model = BDF(debug=False)

        # 11-----10
        # |      |
        # |      |
        # 5------1

        model.add_grid(5, [0., 0., 0.])
        model.add_grid(1, [1., 0., 0.])
        model.add_grid(10, [1., 1., 0.])
        model.add_grid(11, [0., 1., 0.])

        model.add_grid(105, [0.5, 0., 0.])
        model.add_grid(101, [1., 0.5, 0.])
        model.add_grid(110, [0.5, 1., 0.])
        model.add_grid(111, [0., 0.5, 0.])

        #card_lines = ['PLOTEL4', 42,
                      #5, 1, 10, 11]

        card_lines = ['PLOTEL8', 42,
                      5, 1, 10, 11,
                      105, 101, 110, 111]
        model.add_card(card_lines, 'PLOTEL8', comment='plot1')
        model.setup(run_geom_check=True)
        save_load_deck(model)
