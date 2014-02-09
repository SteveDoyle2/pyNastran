import os
import sys
import copy
import shutil
from numpy import allclose

from mapLoads import run_map_loads
from runSpline import run_map_deflections

from pyNastran.utils.log import get_logger
debug = True
log = get_logger(None, 'debug' if debug else 'info')

def run_mapping():
    basepath    = os.path.normpath(os.getcwd())
    configpath  = os.path.join(basepath, 'inputs')
    print "os.getcwd()", os.getcwd()
    globals2 = {}
    locals2 = {}

    fname = "inputs.user_inputs"
    command_module = __import__("inputs.user_inputs", globals2, locals2, ['*'])

    required_inputs = {
        'xref' : None,
        'nastran_call' :  None,
    }
    for key in required_inputs:
        if key not in command_module.__dict__:
            msg = 'fname=%r doesnt contain %r' % (fname, key)
        value = command_module.__dict__[key]
        required_inputs[key] = value
        #print key, value


    nastran_call = required_inputs['nastran_call']
    xref = required_inputs['xref']

    print "nastran_call = %r" % nastran_call
    #print "globals2 =", globals2
    print "xref =", xref
    #print "locals2 =", locals2
    #print "globals()", globals()
    asfd


    workpath    = os.path.join(basepath, 'outputsFinal')

    # load mapping
    cart3dLoads = os.path.join(workpath,  'Cart3d_35000_0.825_10_0_0_0_0.i.triq')
    bdfModel    = os.path.join(configpath,'aeroModel_mod.bdf')
    bdfModelOut = os.path.join(workpath,  'fem_loads_3.bdf')
    # mappingMatrix.new.out - stored in workpath

    # deflection mapping
    cart3dGeom  = os.path.join(configpath, 'Cart3d_bwb.i.tri')
    cart3dGeom2 = os.path.join(workpath, 'Components.i.tri')
    bdf = os.path.join(workpath, 'fem3.bdf')
    op2 = os.path.join(workpath, 'fem3.op2')
    f06 = os.path.join(workpath, 'fem3.f06')

    assert os.path.exists(bdf), '%r doesnt exist' % bdf
    assert os.path.exists(bdfModel), '%r doesnt exist' % bdfModel
    assert os.path.exists(cart3dGeom), '%r doesnt exist' % cart3dGeom
    log.info("basepath = %s" % basepath)

    os.chdir(workpath)
    copyFile(cart3dGeom, 'Components.i.tri')

    nodeList = [20037, 21140, 21787, 21028, 1151, 1886, 2018, 1477, 1023, 1116, 1201, 1116, 1201, 1828, 2589, 1373, 1315, 1571, 1507, 1532, 1317, 1327, 2011, 1445, 2352, 1564, 1878, 1402, 1196, 1234, 1252, 1679, 1926, 1274, 2060, 2365, 21486, 20018, 20890, 20035, 1393, 2350, 1487, 1530, 1698, 1782]
    outfile = open('convergeDeflections.out', 'ab')

    maxADeflectionOld = 0.
    nIterations = 30
    iCart = 1
    for i in range(1, nIterations):
        strI = '_' + str(i)
        assert os.path.exists('Components.i.tri')
        #if i==iCart:
        if 0:
            # run cart3d
            log.info("---running Cart3d #%s---" % i)
            sys.stdout.flush()
            failFlag = os.system('./COMMAND > command.out') # runs cart3d.i.tri, makes Components.i.triq
            assert failFlag == 0, 'Cart3d ./COMMAND failed on iteration #%s' % i
            moveFile('Components.i.triq', cart3dLoads)
            copyFile(cart3dLoads, cart3dLoads+strI)
            copyFile('forces.dat',  'forces.dat'  + strI)
            copyFile('moments.dat', 'moments.dat' + strI)
            copyFile('loadsCC.dat', 'loadsCC.dat' + strI)
            copyFile('history.dat', 'history.dat' + strI)
            os.remove('Components.i.tri') # verifies new Components.i.tri gets created
            sys.stdout.flush()

        # map deflections
        run_map_loads(cart3dLoads, bdfModel, bdfModelOut)  # maps loads
        copyFile(bdfModelOut, bdfModelOut + strI)

        # run nastran
        log.info("---running Nastran #%s---" % i)
        sys.stdout.flush()
        #failFlag = os.system('nastran scr=yes bat=no fem3.bdf') # runs fem3.bdf with fem_loads_3.bdf
        #assert failFlag == 0,'nastran failed on iteration #%s' % i
        copyFile('fem3.op2', 'fem3.op2' + strI)
        copyFile('fem3.f06', 'fem3.f06' + strI)
        os.remove(bdfModelOut) # cleans up fem_loads.bdf

        # map deflections
        (wA, wS) = run_map_deflections(nodeList, bdf, op2, cart3dGeom, cart3dGeom2, log=log)
        assert os.path.exists('Components.i.tri')
        os.remove(op2) # verifies new fem3.op2 was created
        os.remove(f06) # verifies new fem3.f06 was created

        # post-processing
        (maxAID, maxADeflection) = maxDict(wA)
        maxSID = '???'
        maxADeflection = wA[maxAID]
        maxSDeflection = max(wS)[0,0]
        log.info(     "AERO      - i=%s maxAID=%s maxADeflection=%s"   %(i, maxAID, maxADeflection))
        log.info(     "STRUCTURE - i=%s maxSID=%s maxSDeflection=%s"   %(i, maxSID, maxSDeflection))
        outfile.write("AERO      - i=%s maxAID=%s maxADeflection=%s\n" %(i, maxAID, maxADeflection))
        outfile.write("STRUCTURE - i=%s maxSID=%s maxSDeflection=%s\n" %(i, maxSID, maxSDeflection))

        msg  = '\n'+'*'*80+'\n'
        msg += 'finished iteration #%s\n' %(i)
        msg += '*'*80+'\n'
        log.info(msg)

        if allclose(maxADeflection, maxADeflectionOld, atol=0.001):
            break
        maxADeflectionOld = copy.deepcopy(maxADeflection)
        iCart += 1
        sys.stdout.flush()

    outfile.close()
    log.info('---finished runMapping.py---')

def maxDict(dictA):
    k = dictA.keys()
    v = dictA.values()
    max_value = max(v)
    i = v.index(max_value)
    max_key = k[i]
    return (max_key, max_value)

def moveFile(a, b):
    assert a != b, 'a=b=True  a=%r b=%r' % (a, b)
    assert os.path.exists(a),'fileA=%s does not exist...' % a
    if os.path.exists(b):
        os.remove(b)
    shutil.move(a, b)
    assert os.path.exists(b), 'fileB=%r was not moved...' % b

def copyFile(a, b):
    assert a!=b, 'a=b=True  a=%r b=%r' % (a, b)
    assert os.path.exists(a),'fileA=%r does not exist...' % a
    if os.path.exists(b):
        os.remove(b)
    shutil.copyfile(a, b)
    assert os.path.exists(b), 'fileB=%r was not copied...' % b

if __name__=='__main__':
    run_mapping()

