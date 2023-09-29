import unittest
from io import StringIO
import numpy as np

from cpylog import get_logger
from pyNastran.dev.bdf_vectorized3.bdf import BDF, BDFCard, read_bdf#, get_logger2
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
