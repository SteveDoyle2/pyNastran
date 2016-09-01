n = 9 + 1

import struct
from numpy import fromfile

with open('f.bin', 'wb') as f:
    for i in xrange(n):
         f.write(struct.pack('ifi', 1*i, 2.2*i, 3*i))

# this works
dt = [('one', 'int32'),
      ('two', 'float32'),
      ('three', 'int32'),]

with open('f.bin', 'rb') as f:
    #A = fromfile(f, dtype=dt)
    A = fromfile(f, dtype=dt, count=10)
    print(A, len(A)*4)
print("------")


# this fails
if 0:
    with open('f.bin', 'rb') as f:
        d = f.read(n*4)
        print(len(d))
        A = fromfile(d, dtype=dt)
        print(A)
