import numpy as np
from pyNastran.op2.op2 import read_op2
from pyNastran.bdf.utils import parse_patran_syntax_dict

#def test_nodal_avg_stress():
    #op2 = None
    #subcase = None
    #eids = None
    #get_nodal_avg_stress(op2, subcase, eids)
    #get_centroid_stress(op2, subcase, eids)

def get_types_for_eids(bdf, eids):
    elem_prop_types = set()
    for eid in eids:
        elem = bdf.elements[eid]
        etype = elem.type
        ptype = elem.pid.type
        elem_prop_types.add((etype, ptype))
    return elem_prop_types

def get_class_of_elements(etypes):
    if ('CTRIA3', 'PSHELL') in etypes or ('CQUAD4', 'PSHELL') in etypes:
        for etype, ptype in etypes:
            assert etype in ['CTRIA3', 'CQUAD4'], etypes
            assert ptype in ['PSHELL'], etypes
        group = 'shell'
    elif (('CTETRA', 'PSOLID') in etypes or
          ('CPENTA', 'PSOLID') in etypes or
          ('CHEXA', 'PSOLID') in etypes):
        for etype, ptype in etypes:
            assert etype in ['CTETRA', 'CPENTA', 'CHEXA'], etypes
            assert ptype in ['PSOLID'], etypes
        group = 'solid'
    elif ('CBAR', 'PBAR') in etypes or ('CBAR', 'PBARL') in etypes:
        for etype, ptype in etypes:
            assert etype in ['CBAR'], etypes
            assert ptype in ['PBAR', 'PBARL'], etypes
        group = 'bar'
    else:
        raise NotImplementedError(etypes)
    return group

def get_centroid_max_min_principal_stress(bdf, op2, subcase, eids):
    ep_types = get_types_for_eids(bdf, eids)
    group = get_class_of_elements(ep_types)

    eid_max = None
    eid_min = None
    maxp = None
    minp = None
    itime = 0

    if group == 'shell':
        data = []
        imaxp = 6
        iminp = 7

        for obj in [op2.ctria3_stress, op2.cquad4_stress]:
            if subcase in obj:
                headers = obj.get_headers()
                assert headers[imaxp] == 'omax', headers
                assert headers[iminp] == 'omin', headers
                data.append(obj)
        for eid in eids:
            found_eid = False
            for obj in data:
                eids = obj.element_node[:, 0]
                nids = obj.element_node[:, 1]
                assert obj.is_fiber_distance is True, obj.code_information()
                if eid in obj.element:
                    i = np.where(obj.element == eid)[0]
                    print('eid=%s i=%s eid=%s nids=%s' % (eid, i, eids[i], nids[i]))
                    maxi = obj.data[itime, i, imaxp]
                    mini = obj.data[itime, i, iminp]
                    if maxi is None or maxi > maxp:
                        eid_max = eid
                        maxp = maxi
                    if mini is None or mini > minp:
                        eid_min = eid
                        minp = mini
                    found_eid = True
                    break
            assert found_eid is True, eid

    elif group == 'solid':
        data = []
        imaxp = 6
        iminp = 7
        for obj in [op2.ctetra_stress, op2.cpenta_stress, op2.chexa_stress]:
            if subcase in obj:
                headers = obj.get_headers()
                assert headers[imaxp] == 'omax', headers
                assert headers[iminp] == 'omin', headers
                data.append(obj)
                for eid in eids:
                    found_eid = False
                    for obj in data:
                        eids = obj.element_node[:, 0]
                        nids = obj.element_node[:, 1]
                        assert obj.is_fiber_distance is True, obj.code_information()
                        if eid in obj.element:
                            i = np.where(obj.element == eid)[0]
                            #print('eid=%s i=%s eid=%s nids=%s' % (eid, i, eids[i], nids[i]))
                            maxi = obj.data[itime, i, imaxp]
                            mini = obj.data[itime, i, iminp]
                            if maxi is None or maxi > maxp:
                                eid_max = eid
                                maxp = maxi
                            if mini is None or mini > minp:
                                eid_min = eid
                                minp = mini
                            found_eid = True
                            break
                    assert found_eid is True, eid

    #elif group == 'bar':
        #data = []
        #imaxp = 6
        #iminp = 7
        #for obj in [op2.cbar_stress]: # op2.cbar_stress_10nodes
            #if subcase in obj:
                #headers = obj.get_headers()
                #assert headers[imaxp] == 'omax', headers
                #assert headers[iminp] == 'omin', headers
                #data.append(obj)
        #for eid in eids:
            #found_eid = False
            #for obj in data:
                #if eid in obj.element:
                    #found_eid = True
                    #i = np.where(obj.element == eid)[0]
                    #print('eid=%s i=%s' % (eid, i))
            #assert found_eid == True, eid
    else:
        raise NotImplementedError(group)

    return eid_max, maxp, eid_min, minp

#def get_nodal_max_min_principal_stress(bdf, op2, subcase, eids):
    #return eid_max, maxp, eid_min, minp

#def get_nodal_avg_min_principal_stress(bdf, op2, subcase, eids):
    #return eid_max, maxp, eid_min, minp


def main():  # pragma: no cover
    bdf = None
    #op2_filename = 'model.op2'
    op2 = read_op2(op2_filename=None, combine=True, log=None, debug=True,
                   debug_file=None, build_dataframe=False,
                   skip_undefined_matrices=True, mode='msc')

    subcases = [1]

    groups = {
        1 : 'Elm 403082 565514 403084 552195 552196 553965 552204',
    }
    eid_groups = []
    results = {}
    for key, group in sorted(groups.items()):
        eid_group = parse_patran_syntax_dict(group, pound_dict=None)['Elm']
        eid_groups.append(eid_group)
    del groups

    centroid_file = open('centroid.csv', 'w')
    centroid_file.write('group, subcase, eid_max, maxp, eid_min, minp\n')
    for group_id, eids in enumerate(eid_groups):
        # TODO: speed this up by using the same indices
        for subcase in subcases:
            eid_max, maxp, eid_min, minp = get_centroid_max_min_principal_stress(
                bdf, op2, subcase, eids)
            centroid_file.write('%s, %s, %s, %s, %s, %s\n' % (
                group, subcase, eid_max, maxp, eid_min, minp))

    cat('centroid.csv')

def cat(fname):
    with open(fname, 'r') as file_obj:
        print(file_obj.readline().strip())


if __name__ == '__main__':  # pragma: no cover
    main()
