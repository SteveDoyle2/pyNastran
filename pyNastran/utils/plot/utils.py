from typing import Optional
import numpy as np

try:
    import matplotlib.pyplot as plt
    IS_MATPLOTLIB = True
except ModuleNotFoundError:  # pragma: no cover
    IS_MATPLOTLIB = False

def ax2d_to_maximax_minimin(ax: plt.Axes,
                            xlim: Optional[tuple[float, float]]=None,
                            fmt: str='.3g',
                            loc: str='best') -> None:
    xlim = _get_xlim(ax, xlim)
    xmin, xmax = xlim

    ymaxs = []
    ymins = []
    for line in ax.get_lines():
        x = line.get_xdata()
        y = line.get_ydata()
        ix = np.where((xmin <= x) & (x <= xmax))[0]
        #if not len(ix):
            #continue
        yi = y[ix]
        ymaxs.append(yi.max())
        ymins.append(yi.min())

    ymax = max(ymaxs)
    ymin = min(ymins)
    max_label = ('max: %' + fmt) % ymax
    min_label = ('min: %' + fmt) % ymin
    ax.plot(xlim, (ymax, ymax), label=max_label, linestyle='--', color='k')
    ax.plot(xlim, (ymin, ymin), label=min_label, linestyle='--', color='k')
    ax.legend(loc=loc)


def _get_xlim(ax: plt.Axes,
              xlim: Optional[tuple[float, float]] = None) -> tuple[float, float]:
    if xlim is None:
        xmin, xmax = ax.get_xlim()
    else:
        xmin, xmax = xlim

    if xmin is None or xmax is None:
        xmaxs = []
        xmins = []
        for line in ax.get_lines():
            x = line.get_xdata()
            xmaxs.append(x.max())
            xmins.append(x.min())

        # set xmin/xmax to the plot range
        if xmin is None:
            xmin = min(xmins)
        if xmax is None:
            xmax = max(xmaxs)
        xlim = (xmin, xmax)
    return xmin, xmax
