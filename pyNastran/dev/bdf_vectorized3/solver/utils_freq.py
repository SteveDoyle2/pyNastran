import numpy as np
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase


def get_frequencies(model: BDF, subcase: Subcase,
                    omega_ns: np.ndarray) -> np.ndarray:
    """
    Pulls the frequencies from the FREQi cards.

    We need the natural frequencies because some cards
    use that.
    """
    is_freqs = (len(model.frequencies) == 0)
    is_frequency = ('FREQUENCY' not in subcase)
    natural_freq = np.unique(omega_ns) / (2 * np.pi)
    del omega_ns

    if is_freqs or is_frequency:
        nfreq = 1001
        #freq_max = 2 * np.pi * omega_max
        fmax_default = 1.5 * natural_freq.max()

        #freq = np.linspace(0., 1.5*freq_max, num=nfreq)
        frequencies = np.linspace(1., fmax_default, num=nfreq)
    else:
        freq_id, unused_options = subcase['FREQUENCY']
        frequencies_list = []
        for freq in model.frequencies[freq_id]:
            if freq.type in {'FREQ', 'FREQ1', 'FREQ2'}:
                frequencies_list.append(freq.freqs)
            elif freq.type in {'FREQ3', 'FREQ4', 'FREQ5'}:
                freqi = freq.get_frequencies(natural_freq)
                frequencies_list.append(freqi)
            else:  # pragma: no cover
                raise RuntimeError(freq)
        frequencies = np.unique(np.hstack(frequencies_list))
    return frequencies


def slice_freq_set(node_gridtype: np.ndarray,
                   xg: np.ndarray,
                   nnode: int, nfreq: int,
                   node_set: np.ndarray) -> np.ndarray:
    assert xg.ndim == 2, xg.shape
    assert node_gridtype.shape == (nnode, 2), node_gridtype.shape
    assert xg.shape == (nfreq, nnode*6), (xg.shape, (nfreq, nnode*6))
    if node_set[0] != 0:  # 0=all
        assert len(np.unique(node_set)), len(node_set)
        # assert phi.shape == (nnode, nmode), phi.shape
        inode = np.searchsorted(node_gridtype[:, 0], node_set)
        assert np.array_equal(node_gridtype[inode, 0], node_set)
        node_gridtype = node_gridtype[inode, :]
        xg = xg.reshape(nfreq, nnode, 6)[:, inode]
        nnode = len(node_set)
    assert xg.shape == (nfreq, nnode, 6)
    return node_gridtype, xg, nnode
