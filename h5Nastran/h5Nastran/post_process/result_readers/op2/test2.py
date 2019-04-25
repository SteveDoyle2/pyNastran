from pyNastran.op2.dev.op2 import OP2
import struct
import numpy as np


model = OP2(r'example.op2')

print(model._postheaderpos)

directory = model.directory()

for key, data in directory.items():
    print(key, data)


dbnames = model.dbnames

table_ids = dbnames.keys()

dtype = np.dtype([
    ('EKEY', np.int32),
    ('PLY', np.int32),
    ('SX1', np.float32),
    ('SY1', np.float32),
    ('T1', np.float32),
    ('SL1', np.float32),
    ('SL2', np.float32),
    ('A1', np.float32),
    ('MJRP1', np.float32),
    ('MNRP1', np.float32),
    ('TMAX1', np.float32)
])

print(dtype.itemsize)

for table_id in table_ids:
    if table_id != 'OES1C':
        continue
    tables = dbnames[table_id]
    for table in tables:
        fpos, table_len, shape = table
        model._fileh.seek(fpos[0])
        _data = model._fileh.read(864)
        _data = model._fileh.read(524296)
        count = 524296 // 44

        print(count, 44 * count, 524296 - 44 * count)

        data = np.frombuffer(_data, dtype=dtype, count=count)

        _data = model._fileh.read(36)

        print(struct.unpack('9i', _data))

        # print(data)
