from __future__ import print_function
import matplotlib.pyplot as plt
from pyNastran.bdf.bdf import read_bdf
from pyNastran.op2.op2 import OP2

basename = 'transient_test'
#basename = 'test'
model = read_bdf(basename + '.bdf')
tstep = model.tsteps[123]
ndt = sum(tstep.N) + 1

op2 = OP2(debug_file='debug.out')
op2.set_subcases([1])
op2.read_op2(basename + '.op2')

fig, ax = plt.subplots()
ax.plot(op2.accelerations[1].dts)
ax.set_title('index is repeated')
ax.set_xlabel('index location')
ax.set_ylabel('index value')
fig.savefig('index.png')

fig, ax = plt.subplots()
ax.plot(op2.accelerations[1].data[:,0])
ax.set_title('results are duplicated')
fig.savefig('values.png')

print('the length of the vectors should be {}. they are {}'.format(
    ndt,
    len(op2.accelerations[1].dts)),
)

