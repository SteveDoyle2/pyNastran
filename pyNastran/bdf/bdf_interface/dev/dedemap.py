import os
import sys
from collections import defaultdict
bdf_filename = sys.argv[1]
base, ext = os.path.splitext(bdf_filename)
bdf_filename_out = base + '.dedmap' + ext

with open(bdf_filename, 'r') as bdf_file:
    lines = bdf_file.readlines()
nlines = len(lines)

blocks = defaultdict(list)
i = 0
line = lines[i]
while 'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N' not in line:
    print(i, line.strip())
    i += 1
    line = lines[i]

while i < nlines and 'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N' in line:
    line = lines[i]
    print('*', i, line)
    source, dmap = line.strip().split('SUBDMAP  =')
    #dmap = dmap.strip()
    print(f'dmap = {dmap!r}')

    blocks[dmap].append(line)
    i += 1
    line = lines[i]

    blocks[dmap].append(line)
    while 'PAGE' not in line:
        blocks[dmap].append(line)
        i += 1
        line = lines[i]
    for j in range(3):
        i += 1
        line = lines[i]
        print('**', i, line.strip())
    # i += 2
    # line = lines[i]
    # print('**', i, line)
    print('--------------------------')


with open(bdf_filename_out, 'w') as bdf_file_out:
    for block_name, block in blocks.items():
        bdf_file_out.write(f'0        N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N             SUBDMAP  = {block_name}\n')
        bdf_file_out.write(f'      OLD NO.  NEW NO.   ( *I* = INSERTED,   *D* = DELETED )\n')
        for line in block:
            if (line.startswith('     DMAP-DMAP INSTRUCTION') or 
                'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N' in line or 
                'OLD NO.  NEW NO.   ( *I* = INSERTED,   *D* = DELETED )' in line):
                continue
            bdf_file_out.write(line)
print(f'finished writing {bdf_filename_out}')
