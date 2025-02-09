"""
100: neutral file header
   - title
   - version
402: properties
403: Nodes
  - nid, color, pid, type, topology, layer, orient_id, material_orientation,
    geom(5.0), formulation(6.0), contact(7.0), formualtion2(8.0), material_csy, formualtion3(12.0)
404: Elements
413: Layer Data
450: Output Sets
451: Output Data Vectors (obsolete)
601: materials
"""
import os
from collections import defaultdict
from typing import Any
import numpy as np
from cpylog import SimpleLogger
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from pyNastran.op2.tables.oug.oug_displacements import (
    RealDisplacementArray, ComplexDisplacementArray)
# from pyNastran.op2.op2_interface.hdf5_interface import (
#     RealDisplacementArray, ComplexDisplacementArray)
from pyNastran.utils import print_bad_path


block_name_dict = {
    # key, (name, type, read)
    '100': ('Neutral File Header', 'on', True),
    '402': ('Properties', 'other', False),
    '403': ('Nodes', 'geom', True),
    '404': ('Elements', 'geom', True),
    '413': ('Layer Data', 'other', False),
    '450': ('Output Sets', 'results', True),
    '451': ('Output Data Vectors', 'results', True),
    '601': ('Materials', 'other', False),
}

analysis_type_dict = {
    #id: (analysis_type, analysis_step)
    0: ('Unknown', ''),
    1: ('Static', ''),
    2: ('Modes', 'freq='),
    3: ('Transient', 't='),
    4: ('Frequency Response', 't='),
    5: ('Response Spectrum', ''),
    6: ('Random', ''),
    7: ('Linear Buckling', ''),
    8: ('Design Opt', ''),
    9: ('Explicit', ''),
    10: ('Nonlinear Static', ''),
    11: ('Nonlinear Buckling', ''),
    12: ('Nonlinear Transient', 't='),
    19: ('Comp Fluid Dynamics', ''),
    20: ('Steady State Heat Transfer', ''),
    21: ('Transient Heat', 't='),
    22: ('Advanced Nonlinear Static', ''),
    23: ('Advanced Nonlinear Transient', 't='),
    24: ('Advanced Nonlinear Explicit', ''),
    25: ('Static Aeroelasticity', ''),
    26: ('Aerodynamic Flutter', ''),
}

def get_block_lines(lines: list[str],
                    iline: int) -> dict[str, list[str]]:
    block_lines = []
    minus1_count = 0
    while iline < len(lines):
        line = lines[iline]
        if line.startswith('$COM'):
            iline += 1
            continue

        if minus1_count == 2:
            break

        block_lines.append(line.strip())
        if line.strip() == '-1':
            minus1_count += 1
        else:
            minus1_count = 0
        iline += 1
    return block_lines, iline

def get_string(mystr: str) -> str:
    if mystr == '<NULL>':
        return ''
    return mystr


class ScalarResult:
    def __init__(self,
                 header: str,
                 title: str,
                 names: str,
                 results: np.ndarray,
                 location: str):
        self.header = header
        self.title = title
        self.names = names
        self.results = results
        assert location in {'node', 'centroid'}, location


class Neu:
    def __init__(self, node_id=None, xyz=None,
                 elements=None,
                 result_headers=None, results=None,
                 log=None):
        if xyz is None:
            assert node_id is None
            node_id = np.zeros(0, dtype='int32')
            xyz = np.zeros((0, 3), dtype='float64')
        if elements is None:
            elements = {}
        if results is None:
            assert result_headers is None
            result_headers = {}
            results = {}
        #print(xyz.shape)
        #print(f'elements = {elements}')
        #print(list(results.keys()))
        self.log = log
        self.node_id = node_id
        self.xyz = xyz
        self.elements = elements
        self.result_headers = result_headers
        self.results = results
        self.model = to_nastran(self)
        #print(self.model.get_bdf_stats())
        self.organize_results()

    def combine_results_into_sets(self) -> dict[int, Any]:
        set_id_to_results_default = defaultdict(list)

        for (set_id, vector_id), result in self.results.items():
            ncurrent = len(set_id_to_results_default[set_id])
            assert ncurrent + 1 == vector_id
            set_id_to_results_default[set_id].append(result)

        set_id_to_results = {}
        for set_id, results in set_id_to_results_default.items():
            if len(results) == 1:
                set_id_to_results[set_id] = results[0]
            else:
                for (_title, ids, _values) in results:
                    break
                titles = []
                results_list = []
                for (titlei, idsi, valuesi) in results:
                    assert np.array_equal(ids, idsi)
                    if 'Total Translation' in titlei:
                        continue
                    titles.append(titlei)
                    results_list.append(valuesi)
                results2 = np.column_stack(results_list)
                set_id_to_results[set_id] = (titles, idsi, results2)
        return set_id_to_results

    def split_results_by_complex(self, set_id_to_results: dict[int, Any],
                                 ) -> dict[int, Any]:
        log = self.log
        out_result = {}
        for set_id, result_tuple in set_id_to_results.items():
            (names, ids, result) = result_tuple
            self.log.info(f'names = {names}')
            assert isinstance(names, list), names
            found_name = False
            names2 = []
            for iname, name in enumerate(names):
                if ':' in name:
                    resulti = result[:, iname]
                    sline = name.split(':', 1)
                    word = sline[0]
                    if word.endswith('RE'):
                        found_name = True
                        reals = resulti
                        names2.append(name)
                    elif word.endswith('IM'):
                        found_name = True
                        imags = resulti
                        names2.append(name)
                    elif word in {'MAGN', 'PHAS'}:
                        continue
                    else:  # pragma: no cover
                        raise RuntimeError(word)
            if found_name:
                names = names2
                result = reals + 1j * imags
            log.info(f'final names = {names}')
            out_result[set_id] = (names, ids, result)

            ['CPRE:       2M= 0.9k= 0.1',
             'CPIM:       2M= 0.9k= 0.1',
             #'MAGN:       2M= 0.9k= 0.1',
             #'PHAS:       2M= 0.9k= 0.1',
             ]
        return out_result

    def organize_results(self):
        set_id_to_results = self.combine_results_into_sets()
        results = self.split_results_by_complex(set_id_to_results)
        op2_model = OP2(log=self.log)

        ntimes = 1
        ogs = 0
        count = 0
        superelement_adaptivity_index = ''
        pval_step = ''
        sort_method = 0

        for output_set, (names, ids, results) in results.items():
            isubcase = output_set
            #---------------------
            data_dict = self.result_headers[output_set]
            print(data_dict)

            title = data_dict['title']
            if 'T1 Translation' in names[0]:
                n_ids = len(ids)
                assert results.shape == (n_ids, 3), results.shape
                data = np.zeros((ntimes, n_ids, 6))
                data[0, :, :3] = results
                node_gridtype = np.zeros((n_ids, 2), dtype='int32')
                node_gridtype[:, 0] = ids
                table_name = 'OUGV1'
                class_obj = RealDisplacementArray
            else:
                print(data_dict)

            time = data_dict['dt']
            if data_dict['analysis_type'] == 'Transient':
                analysis_code = 6
                # elif trans_word in ['TIME', 'TIME STEP']:  # TODO check name
                #     analysis_code = 6
                # elif trans_word == 'EIGENVALUE':  # normal modes
                #     analysis_code = 2
                # elif trans_word == 'FREQ':  # TODO check name
                #     analysis_code = 5
                # elif trans_word == 'FREQUENCY':
                #     analysis_code = 5
                # elif trans_word == 'COMPLEX EIGENVALUE':
                #     analysis_code = 9

                # {'title': 'Flutter Mode  1    5.44Hz',
                #  'analysis_type': 'Transient',
                #  'from_program': '0',
                #  'analysis_step': 't=5.442459973167886',
                #  'text': (
                #     'ZAERO FLUTTER MODE AT TIME T= 0.18374, FOR PLTFLUT ID=      30',
                #     'FLT FREQ= 5.44246 HZ, FLT SPEED= 903.617, FLT MODE=       1')}
                times = np.array([time], dtype='float64')
                res = class_obj.add_transient_case(
                    table_name, node_gridtype, data, isubcase, times)
                #key = (isubcase, analysis_code, sort_method, count, superelement_adaptivity_index, pval_step)
                key = (isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index, pval_step)
                op2_model.eigenvectors[key] = res
            else:
                if title == 'ZAERO Unsteady Pressure':
                    analysis_code = 9  # complex eigenvalue
                    {'title': 'ZAERO Unsteady Pressure',
                     'dt': 0.1,
                     'analysis_type': 'Unknown',
                     'from_program': '0',
                     'analysis_step': '0.1', 'text': ()}
                    header = title
                    title = 'Cp'
                    res = ScalarResult(header, title, names, results, 'node')
                    key = (isubcase, analysis_code, sort_method, count, ogs, superelement_adaptivity_index, pval_step)
                    op2_model.op2_results.scalars[key] = res
                else:
                    raise RuntimeError(data_dict)
                #print(analysis_type_dict)

        #title_lines = ('ZAERO Unsteady Pressure', ())

        #result_names = [(21, 1), (21, 2), (21, 3), (21, 4)]
        #result_names = [
            # 'CPRE:       2M= 0.9k= 0.1',
            # 'CPIM:       2M= 0.9k= 0.1',
            # 'MAGN:       2M= 0.9k= 0.1',
            # 'PHAS:       2M= 0.9k= 0.1']
        assert len(op2_model.op2_results.scalars) > 0
        self.results_model = op2_model

        title_lines = (
            'Flutter Mode  1    5.44Hz',
            ('ZAERO FLUTTER MODE AT TIME T= 0.00000, FOR PLTFLUT ID=      30',
             'FLT FREQ= 5.44246 HZ, FLT SPEED= 903.617, FLT MODE=       1'))
        # result_names = [
        #     'ZAERO Total Translation',
        #     'ZAERO T1 Translation',
        #     'ZAERO T2 Translation',
        #     'ZAERO T3 Translation']

def to_nastran(model: Neu):
    nastran_model = BDF()
    for nid, xyz in zip(model.node_id, model.xyz):
        nastran_model.add_grid(nid, xyz)

    shell_pids = set()
    for etype, cards in model.elements.items():
        if etype == 'CQUAD4':
            for (eid, pid, nids) in cards:
                nastran_model.add_cquad4(eid, pid, nids)
                shell_pids.add(pid)
        elif etype == 'CTRIA3':
            for (eid, pid, nids) in cards:
                nastran_model.add_ctria3(eid, pid, nids)
                shell_pids.add(pid)
        else:  # pragma: no cover
            raise RuntimeError(etype)

    t = 0.1
    E = 3.0e7
    G = None
    nu = 0.3
    for pid in shell_pids:
        mid = pid
        nastran_model.add_pshell(pid, mid, t)
        nastran_model.add_mat1(mid, E, G, nu)
    return nastran_model


def read_neu(neu_filename: str,
             is_geom: bool=True,
             is_results: bool=True,
             debug: None | bool=True,
             log=None) -> Neu:
    assert os.path.exists(neu_filename), print_bad_path(neu_filename)
    if log is None:
        level = 'warning' if debug is None else (
            'debug' if debug else 'info')
        log = SimpleLogger(level=level)

    with open(neu_filename, 'r') as neu_file:
        lines = neu_file.readlines()

    skip_geom = not is_geom
    skip_results = not is_results
    version, block_dict = read_block_dict(
        log, lines,
        skip_geom=skip_geom,
        skip_results=skip_results,
    )
    node_id = None
    xyz = None
    result_headers = {}
    results = {}
    elements = {}
    #elements = defaultdict(list)
    for key, block in block_dict.items():
        name, result_type, read_flag = block_name_dict.get(key, '???')

        if not read_flag:
            log.warning(f'skipping block {key} ({name})')
            continue

        #print(f'---{key}---')
        #for line in block[:5]:
        #    print(line.rstrip())
        #print('------------')
        if key == '404':
            log.debug(f'reading block {key} ({name})')
            assert len(elements) == 0, elements
            for blocki in block:
                read_elements(elements, blocki, version, debug=False)
            del blocki
        elif key == '403':  # nodes
            log.debug(f'reading block {key} ({name})')
            assert xyz is None, xyz
            nid_list = []
            xyz_list = []
            for blocki in block:
                nidsi, xyzi = read_nodes(blocki, version, debug=False)
                nid_list.append(nidsi)
                xyz_list.append(xyzi)
            node_id = np.hstack(nid_list)
            xyz = np.vstack(xyz_list)
            del nidsi, nid_list, xyzi, xyz_list, blocki
        # elif key == '402':  # properties
        #     pass
        # elif key == '601':  # materials
        #     pass
        # elif key == '413':  # Layer Data
        #     pass
        # elif key == '450':  # Output Sets
        #     pass
        elif key == '450':  # Output Data Name
            log.debug(f'reading block {key} ({name})')
            for iblock, blocki in enumerate(block):
                read_output_headers(
                    log, result_headers,
                    blocki, version, debug=False)
                #print('-------------------')
        elif key == '451':  # Output Data Vectors (obsolete)
            log.debug(f'reading block {key} ({name})')
            assert len(results) == 0, results
            debug = False
            for iblock, blocki in enumerate(block):
                results = read_results(
                    result_headers,
                    blocki, version, debug=debug)
                if debug:
                    print('-------------------')
        else:  # pragma: no cover
            raise RuntimeError(key)

    model = Neu(node_id=node_id, xyz=xyz,
                elements=elements,
                result_headers=result_headers, results=results,
                log=log)
    return model


def read_block_dict(log: SimpleLogger,
                    lines: list[str],
                    skip_geom: bool=False,
                    skip_results: bool=False,
                    ) -> tuple[str, dict]:
    nlines = len(lines)

    minus1_count = 0
    #  -1
    #  100
    #<NULL>
    #12.,
    iline = 0
    # make sure we start at -1
    while iline < nlines:
        line = lines[iline].strip()
        if line == '-1':
            break
        iline += 1
    assert lines[iline].strip() == '-1', lines[0]
    iline += 1

    version = ''
    block_dict = {}
    while iline < nlines:
        line = lines[iline]
        key = line.strip()
        #print(f'{iline}/{nlines}: {line!r} key={key}')
        name, result_type, read_flag = block_name_dict[key]
        if key == '100':
            assert version == '', version
            assert lines[iline+1].strip() == '<NULL>', lines[2]
            version = lines[iline+2].strip(', \n')
            assert version in ['12.', '7.'], f'version={version!r}'
            assert lines[iline+3].strip() == '-1', lines[4]
            assert lines[iline+4].strip() == '-1', lines[5]
            iline += 5
            continue
        elif key in block_name_dict:
            block_lines, iline = get_block_lines(lines, iline)
            if result_type == 'geom' and skip_geom:
                continue
            elif result_type == 'results' and skip_results:
                continue

            if not read_flag:
                log.warning(f'skipping block {key} ({name})')
                continue
            # print(f'-> iline={iline} line={lines[iline]}')
            #print(f'-> iline={iline}')
            #print(block_lines[:3])
            if key not in block_dict:
                block_dict[key] = []
            block_dict[key].append(block_lines)
            continue
        else:  # pragma: no cover
            raise RuntimeError(key)
        iline += 1
    #log.debug(f'made blocks; n={len(block_dict)}')
    return version, block_dict

def int_line(line: str) -> int:
    sline = line.strip(' ,\n').split(',')
    assert len(sline) == 1, sline
    return int(sline[0])

def int_sline(line: str) -> int:
    sline = line.strip(' ,\n').split(',')
    return [int(val) for val in sline]

def float_line(line: str) -> float:
    sline = line.strip(' ,\n').split(',')
    assert len(sline) == 1, sline
    assert 'D' in sline[0]
    return float(sline[0].replace('D', 'E'))

def split_line(line: str) -> list[str]:
    sline = line.strip(' ,\n').split(',')
    return [val.strip() for val in sline]

def read_output_headers(log: SimpleLogger,
                        result_headers: dict[int, dict],
                        block: list[str], version: str,
                        debug: bool=False):
    """
       450
       1,
    Flutter Mode  1    5.44Hz
       0,   3
       5.4424599731678862D+00,
       2,
    ZAERO FLUTTER MODE AT TIME T= 0.00000, FOR PLTFLUT ID=      30
    FLT FREQ= 5.44246 HZ, FLT SPEED= 903.617, FLT MODE=       1
       -1
    """
    assert block[0] == '450', block[0]
    #print(block[0])  # 450
    output_set_str =block[1].strip(' ,\n')
    output_set_id = int(output_set_str)
    #print(f'output_set_id = {output_set_id}')

    # Flutter Mode  1    5.44Hz
    title = block[2].strip()

    from_program, analysis_type_tag = split_line(block[3])
    analysis_type, analysis_step = analysis_type_dict[int(analysis_type_tag)]
    from_program = '???' if from_program == 0 else from_program
    freq = float_line(block[4])

    # if analysis_type == 'Transient' in analysis_type:
    #     log.info(f'program={from_program} analysis_type={analysis_type} time={freq:g}')
    # else:
    #     log.info(f'program={from_program} analysis_type={analysis_type!r}')
    ntext_lines = int_line(block[5])
    #print(f'ntext_lines = {ntext_lines}')
    expected_lines = 8 + ntext_lines
    assert len(block) == expected_lines, (expected_lines, len(block))

    text_lines = ['' if linei == '<NULL>' else linei.strip()
                  for linei in block[6:6+ntext_lines] ]
    text_lines = tuple([line for line in text_lines if line])

    data_dict = {
        'title': title,
        'dt': freq,
        'analysis_type': analysis_type,
        'from_program': from_program,
        'analysis_step': f'{analysis_step}{freq}',
        'text': text_lines,
    }
    #print(f'text_lines = {text_lines}')
    # ZAERO FLUTTER MODE AT TIME T= 0.14699, FOR PLTFLUT ID=      30
    # FLT FREQ= 5.44246 HZ, FLT SPEED= 903.617, FLT MODE=       1
    result_headers[output_set_id] = data_dict

def read_results(headers: dict[int, tuple[str, str, str]],
                 block: list[str], version: str,
                 debug: bool=True):
    """
       1,   1,   1,
    CPRE:       2M= 0.9k= 0.1
          0.,   0.,   0.,
          0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
          0,  0,  3,  8,
          0,  1,  1,
          1001,   1.5713624074338103D-08,

    """
    sline = int_sline(block[1])
    output_set_id, vector_id, one = sline
    assert one == 1, sline
    #print('b:', sline)

    title = block[2]
    #print(f'result title = {title!r}')
    nums = block[8:]
    # print(title)
    ids = []
    values = []
    title_lines = headers[output_set_id]
    if debug:
        print('title_lines =', title_lines)
    #title_end = title_lines[-1]
    results = {
        #(title_end, title): (ids, values),
        (output_set_id, vector_id): (title, ids, values),
    }
    inum = 0
    while inum < len(nums):
        num = nums[inum]
        out = num.strip(', \n').split(',')
        nid_str, *value = [val.strip() for val in out]
        nid = int(nid_str)
        # if nid < 105 or nid > 8503:
        #     print(f'nid = {nid}; values={value}')
        if nid_str == '-1':
            """
                    -1,   0.0000000000000000D+00,
                       1,   2,   1,
                    CPIM:       2M= 0.9k= 0.1
            """
            stop_or_header_line = nums[inum+1]
            inum += 2

            ids = []
            values = []
            if inum >= len(nums):
                break
            #print(f'stop_or_header_line = {stop_or_header_line}')
            sline = int_sline(stop_or_header_line)
            if sline == [-1,]:
                assert inum + 1 == len(nums)
                #print(f'breaking with inum={inum}/{len(nums)}')
                break

            output_set_id, vector_id, one = sline
            assert one == 1, sline

            title = nums[inum]
            #print(f'result title = {title!r}')
            results[(output_set_id, vector_id)] = (title, ids, values)
            # skip over the header
            inum += 6
            continue
        inum += 1

        ids.append(nid_str)
        values.append(value[0].replace('D', 'E'))
        assert len(value) == 1, value
    #assert len(ids) == 0, ids
    #print(len(ids))
    #print(nums[0])
    #print('done with results')
    names = list(results.keys())
    results2 = {}
    for tag, (title, ids, values) in results.items():
        assert len(ids) > 0, ids
        assert len(ids) == len(values)
        ids2 = np.array(ids, dtype='int32')
        values2 = np.array(values, dtype='float64')
        results2[tag] = (title, ids2, values2)
    if debug:
        print(f'result_names = {names}')
    return results2


def read_elements(elements: dict[str, list],
                  block: list[str],
                  version: str, debug: bool=False):
    if version == '7.':
        """
        top=0: line2
        top=2: tri3
        top=3: tri6
        top=4: quad4
        top=5: quad8
        top=6: tet4
        {eid},     {color},     {pid},    {type}, {topology},    {layer},{orient_id},
           1001,       124,         1,        17,         4,         1,         0,         0,         0,         0,         0,         0,
           1002,      1007,      1006,      1001,         0,         0,         0,         0,         0,         0,
              0,         0,         0,         0,         0,         0,         0,         0,         0,         0,
        v:       0.0000000000000000D+00,   0.0000000000000000D+00,   0.0000000000000000D+00,
        offsetA: 0.0000000000000000D+00,   0.0000000000000000D+00,   0.0000000000000000D+00,
        offsetB: 0.0000000000000000D+00,   0.0000000000000000D+00,   0.0000000000000000D+00,
        """
        # for
        elementsi = {}
        cquad4 = []
        ctria3 = []
        ctetra = []
        #cpenta = []
        #chexa = []
        iblock_line = 1
        while iblock_line < len(block) - 2:
            sline1 = block[iblock_line].strip(', \n').split(',')
            sline2 = block[iblock_line+1].strip(', \n').split(',')
            sline3 = block[iblock_line+2].strip(', \n').split(',')
            eid, color, pid, etype, top, *end = sline1
            eid = int(eid)
            pid = int(pid)
            etype = int(etype)
            top = int(top)
            nids_temp = sline2 + sline3
            nids = [int(nid.strip()) for nid in nids_temp[:4]]
            if top == 3:
                top = 'CTRIA3'
                assert len(nids) == 3
                ctria3.append([eid, pid, nids])
            if top == 4:
                top = 'CQUAD4'
                assert len(nids) == 4
                cquad4.append([eid, pid, nids])
            else:  # pragma: no cover
                raise RuntimeError(top)
            if debug:
                print(f'{iblock_line}: eid={eid} pid={pid} etype={etype} top={top}; nids={nids}')
            iblock_line += 7
    else:
        raise RuntimeError(version)
    if cquad4:
        elementsi['CQUAD4'] = cquad4
    if ctria3:
        elementsi['CTRIA3'] = ctria3
    if ctetra:
        elementsi['CTETRA'] = ctetra

    for etype, element_cards in elementsi.items():
        if etype not in elements:
            elements[etype] = []
        elements[etype].extend(element_cards)

def read_nodes(block: list[str], version: str,
               debug: bool=False) -> tuple[np.ndarray, np.ndarray]:
    nids = []
    xyzs = []
    if version == '7.':
        for line in block[1:-2]:
            sline = line.strip('\n ,').split(',')
            # print(f'sline = {sline}')
            # nid,    a,  b,  c, d   e  f   g     h     i     j    x                         y                         z                               k
            # 1001,   0,  0,  1, 46, 0, 0,  0,    0,    0,    0,   8.0000000000000000D+00,   0.0000000000000000D+00,  -5.0000000000000000D+00,         0,
            nid, a, b, c, d, e, f, g, h, i, j, x, y, z, k = sline
            nid = int(nid)
            xyzi = [val.strip().replace('D', 'E') for val in [x, y, z]]
            xyz = np.array(xyzi, dtype='float64')
            if debug:
                print(f'nid={nid} xyz={xyz}')
            nids.append(nid)
            xyzs.append(xyzi)
    else:
        raise RuntimeError(version)
    nid_out = np.array(nids, dtype='int32')
    xyz_out = np.array(xyzs, dtype='float64')
    assert len(nid_out) == len(xyz_out)
    return nid_out, xyz_out
