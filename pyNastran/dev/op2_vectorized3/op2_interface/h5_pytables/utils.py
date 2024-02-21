from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    import pandas as pd

def break_domain_by_case(domains_df: pd.DataFrame, INDEX_DOMAIN) -> pd.DataFrame:
    # these parameters are grouped
    #eigr = mycase['TIME_FREQ_EIGR']
    #eigi = mycase['EIGI']
    #mode = mycase['MODE']

    #idi, subcase, step, analysis, time_freq_eigr, eigi, mode, design_cycle,
    #random, se, afpm, trmc, instance, module, substep, impfid
    #subcase = domains_df[:, 1] # ['SUBCASE']
    subcase = domains_df['SUBCASE']
    isubcase = (subcase != 0)

    # find all unique combinations of these keys
    keys = ['SUBCASE', 'ANALYSIS', 'STEP', 'DESIGN_CYCLE', 'RANDOM', 'SE', 'AFPM', 'TRMC', 'INSTANCE', 'MODULE']
    domain_no_zero = domains_df.loc[isubcase]
    grouped = domain_no_zero.groupby(keys)
    #for g in grouped:
        #print(g)
    grouped_df = grouped.size().reset_index().rename(columns={0:'count'})[keys]
    #step = mycase['STEP']
    #design_cycle = mycase['DESIGN_CYCLE']
    #random = mycase['RANDOM']
    #se = mycase['SE']
    #afpm = mycase['AFPM']
    #trmc = mycase['TRMC']
    #inst = mycase['INSTANCE']
    #module = mycase['MODULE']
    return grouped # _df

def get_name(basename: str, is_freq: bool, subcasei: int, analysis_codei: int,
             modei: int, eigri: float, eigii: float) -> str:
    if modei == 0: # static
        if is_freq:
            name = f'{basename}: subcase={subcasei:d}; freq={eigri:g}; eigi={eigii:g}'
        else:
            name = f'{basename}: subcase={subcasei:d}'
            raise NotImplementedError(analysis_codei)
    elif analysis_codei == 9:
        # we left off eigr/eigi...(eigr=0; eigi!=0)
        name = f'{basename}: subcase={subcasei:d}; loadstep? mode={modei:d} eigi={eigii:g}'
        #raise NotImplementedError(analysis_codei)
    else:
        raise NotImplementedError(analysis_codei)
        #name = f'{basename}: subcase={subcasei:d}; mode={modei:d}; freq={eigri}'
    #name = f'Eigenvector_subcase={subcasei:d}; mode={modei:d}; freq={eigri:g} eigi={eigii:}'
    return name

def get_analysis_code_times(analysis_code, mode, eigr, eigi, idomain):
    is_static = False
    is_modes = False
    is_freq = False
    is_post_buckling = False
    is_transient = False
    is_complex_modes = False

    static_data = None
    modes_data = None
    freq_data = None
    buckling_data = None
    transient_data = None
    complex_mode_data = None
    if analysis_code == 1: # static
        assert mode.max() == 0. and mode.min() == 0., mode
        assert eigr.max() == 0. and eigr.min() == 0., eigr
        assert eigi.max() == 0. and eigi.min() == 0., eigi

        assert len(idomain) == 1, idomain
        static_data = eigr.values
        is_static = True
    elif analysis_code == 2: # modes
        assert eigi.max() == 0. and eigi.min() == 0., eigi
        modes = mode.values # [idomain].values
        eigenvalues = eigr.values # loc[idomain].values
        freq = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
        modes_data = (modes, eigenvalues, freq)
        is_modes = True
    elif analysis_code == 5: # frequency
        assert mode.max() == 0. and mode.min() == 0., mode
        assert eigi.max() == 0. and eigi.min() == 0., eigi
        freq = eigr.loc[idomain].values
        freq_data = freq
        is_freq = True
    elif analysis_code in {6, 1006}: # transient
        assert mode.max() == 0. and mode.min() == 0., mode
        assert eigi.max() == 0. and eigi.min() == 0., eigi
        time = eigr.loc[idomain].values
        transient_data = time
        is_transient = True
    #elif analysis_code == 7: # pre-buckling
    elif analysis_code == 8: # post-buckling
        assert eigi.max() == 0. and eigi.min() == 0., eigi
        modes = mode.values # [idomain].values
        eigenvalues = eigr.values # loc[idomain].values
        freq = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
        buckling_data = (modes, eigenvalues, freq)
        is_post_buckling = True
    elif analysis_code == 9:  # complex eigenvalues
        #assert mode.max() == 0. and mode.min() == 0., mode
        assert eigr.max() == 0. and eigr.min() == 0., eigr
        assert eigi.max() == 0. and eigi.min() == 0., eigi
        complex_mode_data = (mode.values, eigr.values, eigi.values)
        is_complex_modes = True
    else:
        raise NotImplementedError(analysis_code)
    out = (is_static, is_modes, is_freq, is_post_buckling, is_transient, is_complex_modes,
           static_data, modes_data, freq_data, buckling_data, transient_data, complex_mode_data)
    return out
