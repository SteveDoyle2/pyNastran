import h5py
from numpy import linspace, zeros, sin
import matplotlib.pyplot as plt

class Displacement(object):
    def __init__(self, f=None):
        self.f = f

    def build(self, f):
        t = linspace(0., 6., 1000)
        y = sin(t)
        time = f.create_dataset('time', data=t)
        resp = f.create_dataset('response', data=y)
        val = f.create_dataset('value', data=10.0)
        f.value2 = 10.0

def load():
    fname = 'test.hdf5'
    f = h5py.File(fname, 'w')
    d = Displacement()

    f1 = f.create_group('disp1')
    print(f1.name)
    d.build(f1)

    f2 = f.create_group('disp2')
    d.build(f2)
    f.close()
    return fname

def op2():
    hdf5_filename = 'model.hdf5'
    op2_filename = 'model.op2'
    model = OP2()
    model.set_hdf5(hdf5_filename)
    f = model.read_op2(op2_filename, close=True)
    disp = f['displacement'][isubcase]

class OP2_Interface(object):
    def __init__(self, f):
        self.f = f
        self.displacements = {}

    def load(self):
        f = self.f
        for key in f.iterkeys():
            if 'disp' in key:
                print("key =", key)
                i = key[4:]
                group = f[key]
                self.displacements[i] = Displacement(group)
        print(self.displacements.keys())

def op2_interface(f):
    model = OP2_Interface(f)
    model.load()

fname2 = load()
f = h5py.File(fname2, 'r')
print(f.keys())
#print(f['disp1']['time'])
#print(f['disp1']['value'])
#val = f['disp1']['value']
print("value2.name =", f['disp1']['value'].name)
print("value2.value =", f['disp1']['value'].value)
#print("value2.name =", f['disp1']['value'].name)
#print("value2.value =", f['disp1']['value2'].value)
print("keys =", f['disp1'].keys())
#print(val.keys())

op2_interface(f)
#plt.plot(f['disp1']['time'], f['disp1']['response'])
#plt.show()

