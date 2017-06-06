from __future__ import print_function
import os
import sys
import copy
import shutil

from six import iteritems
from six.moves import range
from numpy import allclose

from pyNastran.applications.cart3d_nastran_fsi.map_loads import run_map_loads
from pyNastran.applications.cart3d_nastran_fsi.run_spline import run_map_deflections

from pyNastran.utils.log import get_logger
debug = True
log = get_logger2(None, debug=debug)


def validate_inputs(inputs):
    Mach = inputs['Mach']
    pInf = inputs['pInf']
    qInf = inputs['qInf']
    Sref = inputs['Sref']
    Lref = inputs['Lref']
    xref = inputs['xref']
    isubcase = inputs['isubcase']
    aero_format = inputs['aero_format']

    if not isinstance(isubcase, int):
        raise RuntimeError("isubcase=%r is not an integer greater than 1" % isubcase)
    if not isubcase > 0:
        raise RuntimeError("isubcase=%i is not an greater than 1" % isubcase)

    if not aero_format.lower() in ['cart3d']:
        raise RuntimeError("aero_format=%r is invalid; allowed='cart3d'")
    return True


def load_inputs():
    basepath = os.path.normpath(os.getcwd())
    configpath = os.path.join(basepath, 'inputs')
    workpath = os.path.join(basepath, 'outputsFinal')

    log.info("basepath = %s" % basepath)
    print("os.getcwd() = %s" % os.getcwd())
    globals2 = {}
    locals2 = {}

    fname = "inputs.user_inputs"
    command_module = __import__("inputs.user_inputs", globals2, locals2, ['*'])

    required_inputs = {
        'xref' : None,
        'structural_call' :  None,
        'qInf' : None,
    }
    for key, value in iteritems(command_module.__dict__):
        required_inputs[key] = value

    for key in required_inputs:
        if key not in command_module.__dict__:
            msg = 'fname=%r doesnt contain %r' % (fname, key)
        value = command_module.__dict__[key]
        required_inputs[key] = value
        #print key, value
    required_inputs['configpath'] = configpath
    required_inputs['workpath'] = workpath
    validate_inputs(required_inputs)
    return required_inputs


def run_mapping():
    required_inputs = load_inputs()
    structural_call = required_inputs['structural_call']
    isubcase = required_inputs['isubcase']

    configpath = required_inputs['configpath']
    workpath = required_inputs['workpath']

    print("structural_call = %r" % structural_call)

    # load mapping
    cart3dLoads = os.path.join(workpath, 'Cart3d_35000_0.825_10_0_0_0_0.i.triq')
    bdfModel = os.path.join(configpath, 'aeroModel_mod.bdf')
    bdfModelOut = os.path.join(workpath, 'fem_loads_3.bdf')
    # mappingMatrix.new.out - stored in workpath

    # deflection mapping
    cart3dGeom = os.path.join(configpath, 'Cart3d_bwb.i.tri')
    cart3dGeom2 = os.path.join(workpath, 'Components.i.tri')
    bdf = os.path.join(workpath, 'fem3.bdf')
    #op2 = os.path.join(workpath, 'fem3.op2')
    f06 = os.path.join(workpath, 'fem3.f06')

    assert os.path.exists(bdf), '%r doesnt exist' % bdf
    assert os.path.exists(bdfModel), '%r doesnt exist' % bdfModel
    assert os.path.exists(cart3dGeom), '%r doesnt exist' % cart3dGeom

    os.chdir(workpath)
    copy_file(cart3dGeom, 'Components.i.tri')

    node_list = [
        20037, 21140, 21787, 21028, 1151, 1886, 2018, 1477, 1023, 1116, 1201,
        1116, 1201, 1828, 2589, 1373, 1315, 1571, 1507, 1532, 1317, 1327, 2011,
        1445, 2352, 1564, 1878, 1402, 1196, 1234, 1252, 1679, 1926, 1274, 2060,
        2365, 21486, 20018, 20890, 20035, 1393, 2350, 1487, 1530, 1698, 1782
    ]
    with open('convergeDeflections.out', 'ab') as outfile:
        max_aero_deflection_old = 0.
        niterations = 30
        #icart = 1
        for i in range(1, niterations):
            strI = '_' + str(i)
            assert os.path.exists('Components.i.tri')
            #if i==iCart:
            if 0:
                # run cart3d
                log.info("---running Cart3d #%s---" % i)
                sys.stdout.flush()

                # runs cart3d.i.tri, makes Components.i.triq
                fail_flag = os.system('./COMMAND > command.out')
                assert fail_flag == 0, 'Cart3d ./COMMAND failed on iteration #%s' % i
                move_file('Components.i.triq', cart3dLoads)
                copy_file(cart3dLoads, cart3dLoads + strI)
                copy_file('forces.dat', 'forces.dat' + strI)
                copy_file('moments.dat', 'moments.dat' + strI)
                copy_file('loadsCC.dat', 'loadsCC.dat' + strI)
                copy_file('history.dat', 'history.dat' + strI)
                os.remove('Components.i.tri') # verifies new Components.i.tri gets created
                sys.stdout.flush()

            # map loads
            run_map_loads(required_inputs, cart3dLoads, bdfModel, bdfModelOut)  # maps loads
            copy_file(bdfModelOut, bdfModelOut + strI)

            # run nastran
            log.info("---running Nastran #%s---" % i)
            sys.stdout.flush()
            # runs fem3.bdf with fem_loads_3.bdf
            #fail_flag = os.system('nastran scr=yes bat=no fem3.bdf')
            #assert fail_flag == 0,'nastran failed on iteration #%s' % i
            #copy_file('fem3.op2', 'fem3.op2' + strI)
            copy_file('fem3.f06', 'fem3.f06' + strI)

            # map deflections
            (wA, wS) = run_map_deflections(node_list, bdf, f06, cart3dGeom, cart3dGeom2, log=log)
            #(wA, wS) = run_map_deflections(nodeList, bdf, op2, cart3dGeom, cart3dGeom2, log=log)
            assert os.path.exists('Components.i.tri')

            # cleans up fem_loads.bdf
            os.remove(bdfModelOut)
            #if 0:
                # disabled b/c nastran isn't on this computer
                #os.remove(op2) # verifies new fem3.op2 was created
                #os.remove(f06) # verifies new fem3.f06 was created

            # post-processing
            (max_aero_nid, max_aero_deflection) = max_dict(wA)
            max_structural_nid = '???'
            max_aero_deflection = wA[max_aero_nid]
            max_structural_deflection = max(wS)[0, 0]
            log.info("AERO      - i=%s max_aero_nid=%s max_aero_deflection=%s"   % (
                i, max_aero_nid, max_aero_deflection))
            log.info("STRUCTURE - i=%s max_structural_nid=%s max_structural_deflection=%s"   % (
                i, max_structural_nid, max_structural_deflection))
            outfile.write("AERO      - i=%s max_aero_nid=%s max_aero_deflection=%s\n" % (
                i, max_aero_nid, max_aero_deflection))
            outfile.write("STRUCTURE - i=%s max_structural_nid=%s max_structural_deflection=%s\n" % (
                i, max_structural_nid, max_structural_deflection))

            msg = '\n'+'*' * 80 + '\n'
            msg += 'finished iteration #%s\n' % (i)
            msg += '*' * 80 + '\n'
            log.info(msg)

            if allclose(max_aero_deflection, max_aero_deflection_old, atol=0.001):
                break
            max_aero_deflection_old = copy.deepcopy(max_aero_deflection)
            #icart += 1
            sys.stdout.flush()

    log.info('---finished runMapping.py---')


def max_dict(dictA):
    k = dictA.keys()
    v = dictA.values()
    max_value = max(v)
    i = v.index(max_value)
    max_key = k[i]
    return (max_key, max_value)


def move_file(src, dst):
    assert src != dst, 'a=b=True  src=%r dst=%r' % (src, dst)
    assert os.path.exists(src), 'src=%s does not exist...' % src
    if os.path.exists(dst):
        os.remove(dst)
    shutil.move(src, dst)
    assert os.path.exists(dst), 'fileB=%r was not moved...' % dst


def copy_file(src, dst):
    assert src != dst, 'src=dst=True  src=%r dst=%r' % (src, dst)
    assert os.path.exists(src), 'fileA=%r does not exist...' % src
    if os.path.exists(dst):
        os.remove(dst)
    shutil.copyfile(src, dst)
    assert os.path.exists(dst), 'fileB=%r was not copied...' % dst


if __name__ == '__main__':
    run_mapping()

