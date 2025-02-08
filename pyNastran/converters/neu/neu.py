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
import numpy as np
from pyNastran.utils import print_bad_path
skip_blocks = ['402', '403', '404', '413', '450', '451',
               '601']

def get_block_lines(lines, iline: int):
    block_lines = []
    minus1_count = 0
    while iline < len(lines):
        line = lines[iline]
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

def read_neu(neu_filename: str):
    assert os.path.exists(neu_filename), print_bad_path(neu_filename)
    with open(neu_filename, 'r') as neu_file:
        lines = neu_file.readlines()
    nlines = len(lines)

    minus1_count = 0
    #  -1
    #  100
    #<NULL>
    #12.,
    assert lines[0].strip() == '-1', lines[0]
    iline = 1

    block_dict = {}
    while iline < nlines:
        line = lines[iline]
        key = line.strip()
        #print(f'{iline}/{nlines}: {line!r}')
        if key == '100':
            assert lines[iline+1].strip() == '<NULL>', lines[2]
            version = lines[iline+2].strip(', \n')
            assert version in ['12.', '7.'], f'version={version!r}'
            assert lines[iline+3].strip() == '-1', lines[4]
            assert lines[iline+4].strip() == '-1', lines[5]
            iline += 5
            continue
        elif key in skip_blocks:
            block_lines, iline = get_block_lines(lines, iline)
            # print(f'-> iline={iline} line={lines[iline]}')
            #print(f'-> iline={iline}')
            #print(block_lines[:3])
            assert key not in block_dict, key
            block_dict[key] = block_lines
            continue
        else:
            raise RuntimeError(key)
        iline += 1

    print('made blocks', len(block_dict))
    nodes = []
    for key, block in block_dict.items():
        print(f'---{key}---')
        for line in block[:5]:
            print(line.rstrip())
        print('------------')
        if key == '404':
            read_elements(block, version, debug=False)

        elif key == '403':  # nodes
            xyz = read_nodes(block, version, debug=False)
        elif key == '402':  # properties
            pass
        elif key == '601':  # materials
            pass
        elif key == '413':  # Layer Data
            pass
        elif key == '450':  # Output Sets
            pass
        elif key == '451':  # Output Data Vectors (obsolete)
            read_results(block, version, debug=False)
        else:
            raise RuntimeError(key)

def read_results(block, version, debug=True):
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
    a = block[1]
    title = block[2]
    nums = block[8:]
    # print(a)
    # print(title)
    ids = []
    values = []
    results = {
        title: (ids, values),
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
            inum += 2
            if inum >= len(nums):
                break
            title = nums[inum]
            ids = []
            values = []
            print(f'title = {title!r}')
            results[title] = (ids, values)
            # skip over the header
            inum += 6
            continue
        inum += 1

        ids.append(nid_str)
        values.append(value[0].replace('D', 'E'))
        assert len(value) == 1, value
    print(len(ids))
    print(nums[0])
    print('done with results')
    names = list(results.keys())
    results2 = {}
    for name, (ids, values) in results.items():
        ids2 = np.array(ids, dtype='int32')
        values2 = np.array(values, dtype='float64')
        results2[name] = (ids2, values2)
    print(names)
    return results2


def read_elements(block, version: str, debug: bool=False):
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
        elements = {}
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
            else:
                raise RuntimeError(top)
            print(f'{iblock_line}: eid={eid} pid={pid} etype={etype} top={top}; nids={nids}')
            iblock_line += 7
    else:
        raise RuntimeError(version)
    if cquad4:
        elements['CQUAD4'] = cquad4
    if ctria3:
        elements['CTRIA3'] = ctria3
    if ctetra:
        elements['CTETRA'] = ctetra
    return elements

def read_nodes(block, version: str, debug: bool=False) -> np.ndarray:
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
            xyzs.append(xyzi)
    else:
        raise RuntimeError(version)
    return np.array(xyzs, dtype='float64')

def main():
    neu_filename = r'C:\work\pyNastran\pyNastran\bdf\cards\aero\examples\flutter\case5\cp2anti.neu'
    read_neu(neu_filename)

if __name__ == '__main__':
    main()
