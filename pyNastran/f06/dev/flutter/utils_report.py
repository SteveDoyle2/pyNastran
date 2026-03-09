import os
from pathlib import Path
from collections import defaultdict

from cpylog import SimpleLogger
import numpy as np
import pandas as pd
from pyNastran.utils import PathLike

def write_docx_path(f06_filename: PathLike,
                    ndir_levels: int=1) -> str:
    """add a variable number of directories to the path"""
    nlevels = max(ndir_levels, 0) + 1

    # write the dirname and f06_filename of the file
    parts = Path(f06_filename).parts
    write_parts = parts[-nlevels:]
    path_str = '/'.join(write_parts)
    return path_str


def filenames_to_data_table(filenames2: list[str],
                            ) -> tuple[list[str], list[list[str]]]:
    # remove the extension
    base_filenames = [os.path.splitext(os.path.basename(filename))[0] for filename in filenames2]
    data_table_in = split_by_pattern(base_filenames)
    # print(f'data_table_in = {data_table_in}')
    # find the max length
    lengths = [len(row) for row in data_table_in]
    max_length = max(lengths)
    # print(f'lengths = {lengths}; max={max_length}')

    # define the headers
    # data_row0 = data_table_initial[0]
    headers = [f'{icol + 1:d}'
               for icol in range(max_length)] + ['File']

    # add the filename onto the end
    data_table = [line + ['']*(max_length - len(line)) + [filename]
                  for line, filename in zip(data_table_in, filenames2)]
    return headers, data_table

def split_by_pattern(strings: list[str],
                     delimiter: str='_',
                     group_common: bool=True):
    """
    Main function to split strings based on common pattern.
    """
    if not strings:
        return []

    # Split all strings
    split_lists = [s.split(delimiter) for s in strings]
    if not group_common or len(split_lists) < 2:
        return split_lists
    # split_list0 = split_lists[0]
    # print(f'split_list0 = {split_list0}')

    # Find common prefix
    common_prefix_len = 0
    min_length = min(len(parts) for parts in split_lists)

    for i in range(min_length):
        if all(parts[i] == split_lists[0][i] for parts in split_lists):
            common_prefix_len += 1
            continue
        break
    # print(f'common_prefix_len = {common_prefix_len}')
    # prefix = split_list0[:common_prefix_len]
    # print(f'prefix = {prefix}')

    # Find common suffix
    common_suffix_len = 0
    for i in range(1, min_length - common_prefix_len + 1):
        if all(parts[-i] == split_lists[0][-i] for parts in split_lists):
            common_suffix_len += 1
            continue
        break
    # print(f'common_suffix_len = {common_suffix_len}')

    # Reconstruct
    result = []
    for parts in split_lists:
        new_parts = []
        if common_prefix_len > 0:
            new_parts.append(delimiter.join(parts[:common_prefix_len]))

        middle_start = common_prefix_len
        middle_end = len(parts) - common_suffix_len if common_suffix_len > 0 else len(parts)
        new_parts.extend(parts[middle_start:middle_end])

        if common_suffix_len > 0:
            new_parts.append(delimiter.join(parts[-common_suffix_len:]))
        result.append(new_parts)
    return result


def get_trades(out_table: pd.DataFrame,
               config_text: str,
               log: SimpleLogger) -> tuple[bool, list[str], dict]:
    """
    Returns
    -------
    configs : list[str]
        the first name of the output table
        -> File, but split

    """
    is_passed = True
    columns = out_table.columns
    trade_slines = [val.strip(', ') for val in
                    config_text.strip(';, ').split(';')]

    configs = []
    trades = []
    if 'File' not in columns:
        log.error('Missing "File" from case table')
        is_passed = False
        return is_passed, configs, trades

    nrow = len(out_table)
    files = out_table['File'].to_list()
    if len(trade_slines) == 0 or not is_passed:
        configs = files
        return is_passed, configs, trades

    missing_names = set([])
    for trade in trade_slines:
        if len(trade) == 0:
            continue
        trade_sline = [val.strip() for val in trade.split(',') if val.strip()]
        for col in trade_sline:
            if col not in columns:
                missing_names.add(col)
        icases, icase_dict = get_icases(out_table, trade_sline)
        trades.append((trade_sline, icase_dict))

    if len(missing_names):
        is_passed = False
        missing_names_list = list(missing_names)
        missing_names_list.sort()
        log.error(f'Missing names = {missing_names_list}')

    configs = [''] * nrow
    if not is_passed or len(trades) == 0:
        return is_passed, configs, trades

    cols, icase_dict0 = trades[0]
    depcols = cols[:-1]
    if len(depcols) == 0:
        configs = files
    else:
        config_table = out_table[depcols]
        for i in range(nrow):
            row = config_table.iloc[i, :]
            config = ''
            for col, value in zip(depcols, row):
                config += f'{col}={value}, '
            config = config.strip(', ')
            configs[i] = config

    return is_passed, configs, trades

def get_configs(out_table: pd.DataFrame,
                config_text: str,
                log: SimpleLogger) -> tuple[bool, list[str]]:
    """
    Returns
    -------
    configs : list[str]
        the first name of the output table
        -> File, but split

    """
    is_passed = True
    columns = out_table.columns
    config_keys = config_text.strip(', ').split(',')
    for config_key in config_keys:
        if config_key not in columns:
            log.warning(f'key={config_key!r} is not a column in the case table')

    config_headers_lower = [column.lower().strip() for column in columns]
    config_headers = [column.strip() for column in columns
                      if 'config' in config_headers_lower]

    if 'File' not in columns:
        log.error('missing "File" from case table')
        is_passed = False

    configs = []
    if 'config' in config_headers_lower:
        iconfig_key = config_headers_lower.index('config')
        configs = out_table[iconfig_key].to_list()
        config_headers.remove('config')
    # elif is_passed:
    #     filenames = out_table['File'].tolist()
    #     configs = [os.path.splitext(os.path.basename(filenames))[0]
    #                for filename in filenames]
    else:
        log.error('Missing "Config" from case table')
        is_passed = False

    return is_passed, configs


def data_to_dataframe(headers: list[str],
                      data: list[list[str]],
                      convert_numeric: bool = True) -> pd.DataFrame:
    """Convert numeric columns if requested"""
    df = pd.DataFrame(data, columns=headers)
    if not convert_numeric:
        return df

    for col in df.columns:
        series = df[col]
        words2 = []
        new_words = False
        for word in series.iloc:
            word_lower = word.lower()
            if 'p' in word_lower:
                # 0p2, m0p2
                i = word_lower.index('p')
                word2 = word_lower[:i] + '.' + word_lower[i+1:]
                if word_lower[0] == 'm':
                    word2 = word2[1:]
                # print(f'i={i}; word={word_lower!r} word2={word2!r}')
                try:
                    value2 = float(word2)
                    words2.append(value2)
                    new_words = True
                except ValueError:
                    words2.append(word)
            else:
                words2.append(word)

        # series_floats = series.apply(pd.to_numeric, errors='coerce')
        # if not series_floats.isna().all():
        if new_words:
            df[col] = words2
        else:
            series_floats = series.apply(pd.to_numeric, errors='coerce')
            if not series_floats.isna().all():
                df[col] = series_floats
    return df


def get_icases(df: pd.DataFrame,
               cols: list[str]):
    dfi = df[cols]
    # print(dfi)
    data = dfi.to_numpy()
    nrow = len(data)

    isort = _get_isort_2d_array(data, cols)
    col_data = np.column_stack([np.arange(nrow), data])
    data_sort = col_data[isort, :]
    # print(data_sort)

    case_dict = defaultdict(list)
    for row in data_sort:
        irow, *dep, ind = row
        # print(f'irow={irow} dep={dep} ind={ind}')
        case_dict[tuple(dep)].append(irow)

    icases = []
    for key, icase in case_dict.items():
        # print(f'key={key} icase={icase}')
        icases.append(icase)
    return icases, dict(case_dict)

def _get_isort_2d_array(data: np.ndarray,
                        cols: list[str]) -> np.ndarray:
    """gets the isort flags"""
    nrow = len(data)
    cols_reversed = list(reversed(cols))
    # print(f'cols_reversed = {cols_reversed}')

    isort = np.arange(nrow)
    # print(f'isort = {isort}')
    for icol, col in enumerate(cols_reversed):
        icol_reversed = -icol
        datai = data[:, icol_reversed]
        # print(f'col = {col!r}')
        # print(f'  icol_reversed = {icol_reversed}')
        # print(f'  datai = {datai}')
        isorti = np.argsort(datai)
        # print(f'  isorti = {isorti}')
        isort = isort[isorti]
    # print(f'  isort = {isort}')
    return isort
