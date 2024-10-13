import unittest
import numpy as np
try:
    import matplotlib.pyplot as plt
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

if IS_MATPLOTLIB:
    from pyNastran.utils.plot.utils import ax2d_to_maximax_minimin


class TestPlot(unittest.TestCase):
    def test_plot_time(self):
        fig = plt.figure()
        ax = fig.gca()
        t = np.linspace(0., 1., 100)
        y1 = np.sin(t)
        y2 = np.sin(3*t) * 0.9 - 0.06
        ax.plot(t, y1)
        ax.plot(t, y2)
        ax.grid(True)
        ax2d = ax2d_to_maximax_minimin(ax)

    def test_plot_freq(self):
        fig = plt.figure()
        ax = fig.gca()
        freq = np.logspace(1, 3., num=100)
        Qs = np.array([5., 10., 15.])
        freq_naturals = [20, 100, 200]
        for Q in Qs:
            y = sdof_response(
                freq, freq_naturals=freq_naturals, Q=Q)
            ax.loglog(freq, y)
        ax.grid(True)
        xlim = None
        xlim = (None, 600.)
        #ax.set_xlim(xlim)
        #ax.set_ylim([0.9, 1.2])
        ax2d = ax2d_to_maximax_minimin(
            ax, fmt='.3f', xlim=xlim)

        plt.show()


def sdof_response(freq: np.ndarray,
                  freq_naturals: float | np.ndarray,
                  Q: float=10.0) -> np.ndarray:
    zeta = 1 / (2 * Q)
    if isinstance(freq_naturals, np.ndarray):
        pass
    elif isinstance(freq_naturals, float):
        freq_naturals = np.array([freq_naturals])
    else:
        assert isinstance(freq_naturals, (list, tuple))
        freq_naturals = np.asarray(freq_naturals)
    assert freq_naturals.ndim == 1, freq_naturals.shape

    rho = freq[:, np.newaxis] / freq_naturals[np.newaxis, :]
    tzw = (2 * zeta * rho)**2
    num = 1 + tzw
    denom = (1 - rho**2)**2 + tzw
    y = np.sqrt(num / denom)
    return y
