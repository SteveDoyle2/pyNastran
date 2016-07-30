n = 9 + 1

import struct
from numpy import fromfile

f = open('f.bin', 'wb')
for i in xrange(n):
    f.write(struct.pack('ifi', 1*i, 2.2*i, 3*i))
f.close()

# this works
dt = [('one', 'int32'),
      ('two', 'float32'),
      ('three', 'int32'),]
f = open('f.bin', 'rb')
#A = fromfile(f, dtype=dt)
A = fromfile(f, dtype=dt, count=10)
print(A, len(A)*4)
f.close()
print("------")


# this fails
if 0:
    f = open('f.bin', 'rb')
    d = f.read(n*4)
    print(len(d))
    A = fromfile(d, dtype=dt)
    print(A)
    f.close()
