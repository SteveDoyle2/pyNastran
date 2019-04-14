from collections import defaultdict
from numpy import amax, amin, mean, array, sqrt  #, sum
from numpy.linalg import eigvalsh  # type: ignore

def vonMises3D(o1, o2, o3):
    ovm = sqrt(0.5 * (o1-o2)**2 + (o2-o3)**2 + (o1-o3)**2)
    return ovm

def vonMises2D(o1, o2):
    return sqrt(0.5 * o1**2 + o2**3 - o1*o2)

def get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options=None):
    """
    Supports:
    - plateStress
    - solidStress
    - compositePlateStress (NA)

    options = {
        'mode': 'derive/avg',  # derive/avg, avg/derive
        'layers' : 'max',      # max, min, avg
       #'ilayers' : None,      # None or [1, 2, ... N]
        'location' : 'node',   # node, centroid
    }

    .. todo:  this isn't part of OP2() because it's not done
    .. todo:  doesn't support transient, frequency, real/imaginary data
    .. todo:  add 'sum', 'difference' for 'layers'?
    .. todo:  hasn't been tested

    """
    raise NotImplementedError()
    #assert options['mode'] in ['derive/avg', 'avg/derive'], options['mode']
    #assert options['layers'] in ['max', 'min', 'avg'], options['layers']
    #assert options['location'] in ['node', 'centroid'], options['location']
    #assert options['mode'] in ['derive/avg', 'avg/derive'], options['mode']
    #assert options['mode'] in ['derive/avg', 'avg/derive'], options['mode']

    #layer_map = {
        #'max': amax,
        #'min': amin,
        #'avg': mean,
        ##'sum': sum,
    #}
    #mode = options['mode']
    #layer_func = layer_map[options['layers']]
    #location = options['location']

    #results = {
        #'x': defaultdict(list),
        #'y': defaultdict(list),
        #'z': defaultdict(list),
        #'xy': defaultdict(list),
        #'yz': defaultdict(list),
        #'xz': defaultdict(list),
        #'maxP': defaultdict(list),
        #'minP': defaultdict(list),
        #'vonMises': defaultdict(list),  # 3D von mises
        #'vonMises2D': defaultdict(list),  # 3D von mises
    #}
    #if isubcase in model.solidStress:
        #case = model.solidStress[isubcase]
        #if case.is_von_mises:
            #vmWord = 'vonMises'
        #else:
            #vmWord = 'maxShear'
        #assert vmWord == 'vonMises', vmWord

        #if location == 'node':  # derive/avg
            #for eid in case.ovmShear:
                #node_ids = eid_to_nid_map[eid]
                #for nid in node_ids:
                    #results['x'   ][nid].append(case.oxx[eid][nid])
                    #results['y'   ][nid].append(case.oyy[eid][nid])
                    #results['z'   ][nid].append(case.ozz[eid][nid])
                    #results['xy'  ][nid].append(case.txy[eid][nid])
                    #results['yz'  ][nid].append(case.tyz[eid][nid])
                    #results['xz'  ][nid].append(case.txz[eid][nid])
                    #results['maxP'][nid].append(case.o1[eid][nid])
                    #results['minP'][nid].append(case.o3[eid][nid])
                    #results['vonMises'][nid].append(case.ovmShear[eid][nid])
        #elif location == 'centroid':
            #for eid in case.ovmShear:  # derive/avg
                #node_ids = eid_to_nid_map[eid]
                #for nid in node_ids:
                    #results['x'   ][nid].append(case.oxx[eid]['CENTER'])
                    #results['y'   ][nid].append(case.oyy[eid]['CENTER'])
                    #results['z'   ][nid].append(case.ozz[eid]['CENTER'])
                    #results['xy'  ][nid].append(case.txy[eid]['CENTER'])
                    #results['yz'  ][nid].append(case.tyz[eid]['CENTER'])
                    #results['xz'  ][nid].append(case.txz[eid]['CENTER'])
                    #results['maxP'][nid].append(case.o1[eid]['CENTER'])
                    #results['minP'][nid].append(case.o3[eid]['CENTER'])
                    #results['vonMises'][nid].append(case.ovmShear[eid]['CENTER'])
        #else:
            #raise RuntimeError('location=%r' % location)

    #if isubcase in model.plateStress:
        #case = model.plateStress[isubcase]
        #if case.nonlinear_factor is not None: # transient
            #return
        #if case.is_von_mises:
            #vmWord = 'vonMises'
        #else:
            #vmWord = 'maxShear'

        #assert vmWord == 'vonMises', vmWord
        #if location == 'node':
            #for eid in case.ovmShear:
                #node_ids = eid_to_nid_map[eid]
                #eType = case.eType[eid]
                #if eType in ['CQUAD4', 'CQUAD8']:
                    ##cen = 'CEN/%s' % eType[-1]
                    #assert len(node_ids[:4]) == 4, len(node_ids[:4])
                    #if node_ids[0] in case.oxx[eid]:
                        ## bilinear
                        #for nid in node_ids[:4]:
                            #results['x'   ][nid].append(layer_func(case.oxx[eid][nid]))
                            #results['y'   ][nid].append(layer_func(case.oyy[eid][nid]))
                            #results['xy'  ][nid].append(layer_func(case.txy[eid][nid]))
                            #results['maxP'][nid].append(layer_func(case.majorP[eid][nid]))
                            #results['minP'][nid].append(layer_func(case.minorP[eid][nid]))
                            #results['vonMises'][nid].append(layer_func(case.ovmShear[eid][nid]))
                    #else:
                        ##cen = 'CEN/%s' % eType[-1]
                        #cen = 0
                        #for nid in node_ids[:4]:
                            #results['x'   ][nid].append(layer_func(case.oxx[eid][cen]))
                            #results['y'   ][nid].append(layer_func(case.oyy[eid][cen]))
                            #results['xy'  ][nid].append(layer_func(case.txy[eid][cen]))
                            #results['maxP'][nid].append(layer_func(case.majorP[eid][cen]))
                            #results['minP'][nid].append(layer_func(case.minorP[eid][cen]))
                            #results['vonMises'][nid].append(layer_func(case.ovmShear[eid][cen]))
                #elif eType in ['CTRIA3', 'CTRIA6']:
                    ##cen = 'CEN/%s' % eType[-1]
                    #cen = 0
                    #assert len(node_ids[:3]) == 3, len(node_ids[:3])
                    #for nid in node_ids[:3]:
                        #results['x'   ][nid].append(layer_func(case.oxx[eid][cen]))
                        #results['y'   ][nid].append(layer_func(case.oyy[eid][cen]))
                        #results['xy'  ][nid].append(layer_func(case.txy[eid][cen]))
                        #results['maxP'][nid].append(layer_func(case.majorP[eid][cen]))
                        #results['minP'][nid].append(layer_func(case.minorP[eid][cen]))
                        #results['vonMises'][nid].append(layer_func(case.ovmShear[eid][cen]))
                #else:
                    #raise NotImplementedError(eType)
        #elif location == 'centroid':
            #for eid in case.ovmShear:
                #node_ids = eid_to_nid_map[eid]
                #eType = case.eType[eid]
                #if eType in ['CQUAD4', 'CQUAD8', 'CTRIA3', 'CTRIA6']:
                    ##cen = 'CEN/%s' % eType[-1]
                    ## cen  0
                    #pass
                #else:
                    #raise NotImplementedError(eType)
                #for nid in node_ids:
                    #results['x'   ][nid].append(layer_func(case.oxx[eid][cen]))
                    #results['y'   ][nid].append(layer_func(case.oyy[eid][cen]))
                    #results['xy'  ][nid].append(layer_func(case.txy[eid][cen]))
                    #results['maxP'][nid].append(layer_func(case.majorP[eid][cen]))
                    #results['minP'][nid].append(layer_func(case.minorP[eid][cen]))
                    #results['vonMises'][nid].append(layer_func(case.ovmShear[eid][cen]))
        #else:
            #raise RuntimeError('location=%r' % location)

    #if isubcase in model.compositePlateStress:
        #case = model.compositePlateStress[isubcase]
        #if case.nonlinear_factor is not None: # transient
            #return
        #if case.is_von_mises:
            #vmWord = 'vonMises'
        #else:
            #vmWord = 'maxShear'

        #assert vmWord == 'vonMises', vmWord
        #if location == 'node':
            #for eid in case.ovmShear:
                #node_ids = eid_to_nid_map[eid]
                #eType = case.eType[eid]
                #if eType in ['CQUAD4', 'CQUAD8']:
                    #assert len(node_ids[:4]) == 4, len(node_ids[:4])
                    #for nid in node_ids[:4]:
                        #results['x'   ][nid].append(layer_func(case.o11[eid]))
                        #results['y'   ][nid].append(layer_func(case.o22[eid]))
                        #results['xy'  ][nid].append(layer_func(case.t12[eid]))
                        #results['maxP'][nid].append(layer_func(case.majorP[eid]))
                        #results['minP'][nid].append(layer_func(case.minorP[eid]))
                        #results['vonMises'][nid].append(layer_func(case.ovmShear[eid]))
                #elif eType in ['CTRIA3', 'CTRIA6']:
                    #cen = 'CEN/%s' % eType[-1]
                    #assert len(node_ids[:3]) == 3, len(node_ids[:3])
                    #for nid in node_ids[:3]:
                        #results['x'   ][nid].append(layer_func(case.o11[eid]))
                        #results['y'   ][nid].append(layer_func(case.o22[eid]))
                        #results['xy'  ][nid].append(layer_func(case.t12[eid]))
                        #results['maxP'][nid].append(layer_func(case.majorP[eid]))
                        #results['minP'][nid].append(layer_func(case.minorP[eid]))
                        #results['vonMises'][nid].append(layer_func(case.ovmShear[eid]))
                #else:
                    #raise NotImplementedError(eType)
        #elif location == 'centroid':
            #for eid in case.ovmShear:
                #node_ids = eid_to_nid_map[eid]
                #eType = case.eType[eid]
                #for nid in node_ids:
                    #results['x'   ][nid].append(layer_func(case.o11[eid]))
                    #results['y'   ][nid].append(layer_func(case.o22[eid]))
                    #results['xy'  ][nid].append(layer_func(case.t12[eid]))
                    #results['maxP'][nid].append(layer_func(case.majorP[eid]))
                    #results['minP'][nid].append(layer_func(case.minorP[eid]))
                    #results['vonMises'][nid].append(layer_func(case.ovmShear[eid]))
        #else:
            #raise RuntimeError('location=%r' % location)

    #if mode == 'derive/avg':
        #for result_name, result in results.items():
            #for nid, datai in result.items():
                #results[result_name][nid] = mean(datai)
    #elif mode == 'avg/derive':
        #for result_name in ['x', 'y', 'z', 'xy', 'yz', 'xz']:
            #for nid, datai in results[result_name].items():
                #results[result_name][nid] = mean(datai)

        #for nid in results['maxP']:
            #oxx = results['x'][nid]
            #oyy = results['y'][nid]
            #ozz = results['z'][nid]

            #txy = results['xy'][nid]
            #tyz = results['yz'][nid]
            #txz = results['xz'][nid]

            #if not isinstance(ozz, float):
                #ozz = 0.
            #if not isinstance(txy, float):
                #txy = 0.
            #if not isinstance(tyz, float):
                #tyz = 0.
            #if not isinstance(txz, float):
                #txz = 0.

            ## 3D
            #A = array([
                #[oxx, txy, txz],
                #[ 0., oyy, tyz],
                #[ 0.,  0., ozz],
            #])
            #eigs = eigvalsh(A, UPLO='U')
            #maxP = eigs.max()
            #minP = eigs.min()
            #results['maxP'][nid] = maxP
            #results['minP'][nid] = minP

            ## 2D
            #A2 = array([
                #[oxx, txy],
                #[ 0., oyy],
            #])
            #eigs2 = eigvalsh(A2, UPLO='U')
            ##maxP2 = eigs2.max()
            ##minP2 = eigs2.min()

            #results['vonMises'][nid] = vonMises3D(*eigs)
            #results['vonMises2D'][nid] = vonMises2D(*eigs2)
    #else:
        #raise RuntimeError('mode=%r' % mode)
    #return results

def solid_bending():
    model_path = os.path.join(pyNastran.__path__[0], '..', 'models')
    bdf_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.bdf')
    op2_filename = os.path.join(model_path, 'solid_bending', 'solid_bending.op2')

    options1 = {
        'location': 'node',
        'mode': 'derive/avg',
        'layers': 'max',
    }
    isubcase = 1
    eid_to_nid_map, model = setup(bdf_filename, op2_filename)
    vm_derive_avg = get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options1)['vonMises']

    options2 = {
        'location': 'node',
        'mode': 'avg/derive',
        'layers': 'max',
    }
    vm_avg_derive = get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options2)['vonMises']
    for nid, stressi in sorted(vm_derive_avg.items()):
        print("nid=%-3s derive/avg=%-8g avg/derive=%-8g diff=%g" % (nid, stressi, vm_avg_derive[nid], stressi - vm_avg_derive[nid]))

def sol_101_elements():
    model_path = os.path.join(pyNastran.__path__[0], '..', 'models')
    bdf_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.bdf')
    op2_filename = os.path.join(model_path, 'sol_101_elements', 'static_solid_shell_bar.op2')

    options1 = {
        'location': 'centroid',
        'mode': 'derive/avg',
        'layers': 'min',
    }
    isubcase = 1
    eid_to_nid_map, model = setup(bdf_filename, op2_filename)
    vm_derive_avg = get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options1)['vonMises']

    options2 = {
        'location': 'centroid',
        'mode': 'avg/derive',
        'layers': 'min',
    }
    vm_avg_derive = get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options2)['vonMises']
    for nid, stressi in sorted(vm_derive_avg.items()):
        print("nid=%-3s derive/avg=%-8g avg/derive=%-8g diff=%g" % (nid, stressi, vm_avg_derive[nid], stressi - vm_avg_derive[nid]))

def setup(bdf_filename, op2_filename):
    fem = BDF(debug=False)
    fem.read_bdf(bdf_filename)

    eid_to_nid_map = {}
    for eid, element in fem.elements.items():
        node_ids = element.node_ids
        eid_to_nid_map[eid] = node_ids
    del fem

    model = OP2()
    model.set_results(['solidStress', 'plateStress'])
    model.read_op2(op2_filename)
    return eid_to_nid_map, model

def run(eid_to_nid_map, model, isubcase, options):
    stress = get_nodal_averaged_stress(model, eid_to_nid_map, isubcase, options)
    vm = stress['vonMises']
    del stress
    return vm

if __name__ == '__main__':  # pragma: no cover
    import os
    import pyNastran
    from pyNastran.bdf.bdf import BDF
    from pyNastran.op2.op2 import OP2
    sol_101_elements()
    #solid_bending()
