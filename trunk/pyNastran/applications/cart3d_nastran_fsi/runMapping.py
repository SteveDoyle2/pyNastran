import os
import sys
import copy
import shutil
from numpy import allclose

from mapLoads import runMapLoads
from runSpline import run as runMapDeflections

from logger import dummyLogger
loggerObj = dummyLogger()
log = loggerObj.startLog('debug') # or info

def run():
    basepath    = os.path.normpath(os.getcwd())
    configpath  = os.path.join(basepath,'inputs')
    workpath    = os.path.join(basepath,'outputsFinal')

    # load mapping
    cart3dLoads = os.path.join(workpath,  'Cart3d_35000_0.825_10_0_0_0_0.i.triq')
    bdfModel    = os.path.join(configpath,'aeroModel_mod.bdf')
    bdfModelOut = os.path.join(workpath,  'fem_loads_3.bdf')
    # mappingMatrix.new.out - stored in workpath

    # deflection mapping
    cart3dGeom  = os.path.join(configpath,'Cart3d_bwb.i.tri')
    cart3dGeom2 = os.path.join(workpath,'Components.i.tri')
    bdf = os.path.join(workpath,'fem3.bdf')
    op2 = os.path.join(workpath,'fem3.op2')
    f06 = os.path.join(workpath,'fem3.f06')

    assert os.path.exists(bdf       ),'|%s| doesnt exist' %(bdf)
    assert os.path.exists(bdfModel  ),'|%s| doesnt exist' %(bdfModel)
    assert os.path.exists(cart3dGeom),'|%s| doesnt exist' %(cart3dGeom)
    log().info("basepath = %s" %(basepath))

    os.chdir(workpath)
    copyFile(cart3dGeom, 'Components.i.tri')

    nodeList = [20037, 21140, 21787, 21028, 1151, 1886, 2018, 1477, 1023, 1116, 1201, 1116, 1201, 1828, 2589, 1373, 1315, 1571, 1507, 1532, 1317, 1327, 2011, 1445, 2352, 1564, 1878, 1402, 1196, 1234, 1252, 1679, 1926, 1274, 2060, 2365, 21486, 20018, 20890, 20035, 1393, 2350, 1487, 1530, 1698, 1782]
    outfile = open('convergeDeflections.out','ab')
    
    maxADeflectionOld = 0.
    nIterations = 30
    iCart = 1
    for i in range(1,nIterations):
        strI = '_'+str(i)
        assert os.path.exists('Components.i.tri')
        #if (i==iCart):
        if (1):
            # run cart3d
            log().info("---running Cart3d #%s---" %(i))
            sys.stdout.flush()
            failFlag = os.system('./COMMAND > command.out') # runs cart3d.i.tri, makes Components.i.triq
            assert failFlag == 0,'cart3d failed on iteration #%s' %(i)
            moveFile('Components.i.triq', cart3dLoads)
            copyFile(cart3dLoads, cart3dLoads+strI)
            copyFile('forces.dat', 'forces.dat' +strI)
            copyFile('moments.dat','moments.dat'+strI)
            copyFile('loadsCC.dat','loadsCC.dat'+strI)
            copyFile('history.dat','history.dat'+strI)
            os.remove('Components.i.tri') # verifies new Components.i.tri gets created
            sys.stdout.flush()

        # map deflections
        runMapLoads(cart3dLoads,bdfModel,bdfModelOut)  # maps loads
        copyFile(bdfModelOut, bdfModelOut+strI)

        # run nasran
        log().info("---running Nastran #%s---" %(i))
        sys.stdout.flush()
        failFlag = os.system('nastran scr=yes bat=no fem3.bdf') # runs fem3.bdf with fem_loads_3.bdf
        assert failFlag == 0,'nastran failed on iteration #%s' %(i)
        copyFile('fem3.op2', 'fem3.op2'+strI)
        copyFile('fem3.f06', 'fem3.f06'+strI)
        os.remove(bdfModelOut) # cleans up fem_loads.bdf
        
        # map deflections
        (wA,wS) = runMapDeflections(nodeList,bdf,f06,cart3dGeom,cart3dGeom2)
        assert os.path.exists('Components.i.tri')
        os.remove(op2) # verifies new fem3.op2 was created
        os.remove(f06) # verifies new fem3.f06 was created
        
        # post-processing
        (maxAID,maxADeflection) = maxDict(wA)
        maxSID = '???'
        maxADeflection = wA[maxAID]
        maxSDeflection = max(wS)[0,0]
        log().info(   "AERO      - i=%s maxAID=%s maxADeflection=%s"   %(i,maxAID,maxADeflection))
        log().info(   "STRUCTURE - i=%s maxSID=%s maxSDeflection=%s"   %(i,maxSID,maxSDeflection))
        outfile.write("AERO      - i=%s maxAID=%s maxADeflection=%s\n" %(i,maxAID,maxADeflection))
        outfile.write("STRUCTURE - i=%s maxSID=%s maxSDeflection=%s\n" %(i,maxSID,maxSDeflection))

        msg  = '\n'+'*'*80+'\n'
        msg += 'finished iteration #%s\n' %(i)
        msg += '*'*80+'\n'
        log().info(msg)

        if allclose(maxADeflection,maxADeflectionOld,atol=0.001):
            break
        maxADeflectionOld = copy.deepcopy(maxADeflection)
        iCart+=1
        sys.stdout.flush()
    ###
    outfile.close()
    log().info('---finished runMapping.py---')

def maxDict(dictA):
    k = dictA.keys()
    v = dictA.values()
    maxValue = max(v)
    
    for key,value in dictA.items():
        if value==maxValue:
            maxKey = key
            break
        ###
    ###
    return (maxKey,maxValue)

def moveFile(a,b):
    assert a!=b, 'a=b=True  a=|%s| b=|%s|' %(a,b)
    assert os.path.exists(a),'fileA=%s does not exist...' %(a)
    if os.path.exists(b):
        os.remove(b)
    shutil.move(a,b)
    assert os.path.exists(b),'fileB=%s was not moved...' %(b)

def copyFile(a,b):
    assert a!=b, 'a=b=True  a=|%s| b=|%s|' %(a,b)
    assert os.path.exists(a),'fileA=%s does not exist...' %(a)
    if os.path.exists(b):
        os.remove(b)
    shutil.copyfile(a,b)
    assert os.path.exists(b),'fileB=%s was not copied...' %(b)

if __name__=='__main__':
    run()

