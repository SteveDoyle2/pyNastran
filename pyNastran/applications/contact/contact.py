from __future__ import division
import os
import copy

from numpy import zeros, ones, array, where, allclose

from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.fieldWriter import print_card
#from pyNastran.bdf.writePath import write_include


max_deflection_error = 0.01

def split_model(model, nids, func=None):
    if func is None:
        def func(element):
            return True
    nidi = max(model.nodes) + 1
    old_nids = set(model.nodes.keys())

    nodes_map = {}
    for nid in nids:
        node = model.nodes[nid]
        node2 = copy.deepcopy(node)
        node2.nid = nidi
        model.nodes[nidi] = node2
        nodes_map[nid] = nidi
        nidi += 1
    
    for eid, element in model.elements.iteritems():
        #print dir(element)
        e_nids = element.nodeIDs()
        e_nids2 = [nid for nid in e_nids if nid in nids]
        #e_nids2 = [nid if nid in nids else False for nid in e_nids]

        if element.type in ['CTRIA3', 'CQUAD4']:
            #print "e_nids2 =", e_nids2
            if func(element):
                for nid in e_nids2:
                    i = e_nids.index(nid)

                    nidi = nodes_map[nid]
                    #print "    nid=%s -> nidi=%s" % (nid, nidi)
                    node2 = model.nodes[nidi]
                    #print node2
                    element.nodes[i] = node2
            else:
                continue
                    
        else:
            raise NotImplementedError(element.type)

    all_nids = set([])
    for eid, element in model.elements.iteritems():
        e_nids = element.nodeIDs()
        for nid in e_nids:
            all_nids.add(nid)
    lost_nids = set(model.nodes.keys()) - all_nids
    #all_nids = list(all_nids)
    
    for nid in lost_nids:
        del model.nodes[nid]
    
    #print "all_nids", all_nids
    #print "old_nids", old_nids
    #new_nids = all_nids.remove(old_nids.intersection(lost_nids) )
    #print "new_nids", new_nids
    #nid0 = min(new_nids)
    #model.nodes[nid]._comment = 'updated...\n'
    

def run(contact_surfaces, main_bdf, main_op2):
    contact_bdf = 'spring.bdf'
    subcase_id = 1

    eid_groups, nodes_groups, stiffnesses, errors = setup_contact(main_bdf, contact_bdf, contact_surfaces)
    run_nastran(main_bdf)
    nerrors = []
    ierrrors = None
    ierrors_old = parse_op2(contact_bdf, op2_name, subcase_id, eid_groups, nodes_groups, stiffnesses, errors)

    while 1:
        error = parse_op2(contact_bdf, op2_name, subcase_id, eid_groups, nodes_groups, stiffnesses, errors)
        errors.append(error)
        if error == 0 or error > error_old:
            break
        error_old = error


def run_nastran(main_bdf):
    os.system("nastran " + main_bdf + " scr=yes")


def update_nodes_cid(model, nids_group, cid):
    nnodes = len(nids_group)
    nodes_left = zeros((nnodes, 3), 'float64')

    for i, nid in enumerate(nids_group):
        node = model.nodes[nid]
        if node.cid != 0:
            xyz = node.PositionWRT(cid=cid)
            node.UpdateCoord(xyz, cid)
        nodes[i, :] = node.Position()
    return nodes

def get_max_eid(model):
    """
    assumes no RBE2/3
    """
    ielement_max = max(model_main.elements)
    return max([ielement_max])

def setup_contact(main_bdf, contact_bdf, contact_surfaces):
    main_bdf2 = 'main2.bdf'

    model_main = BDF()
    model_main.read_bdf(main_bdf)

    # ------------------------------------------------------------------------
    # update the case control deck

    cc = model.caseControlDeck
    subcase_ids = cc.get_subcase_list()
    subcase0_id = subcase_ids[0]
    #subcases = model.caseControlDeck.get_subcases()
    #assert len(subcases) == 1, len(subcases)
    #subcase0_id = subcases.keys()[0]
    #subcase = subcases[subcase0_id]
    #subcase['FORCE'] = ['FORCE', 'ALL']
    
    #cc.get_subcase_parameter(subcase0_id, 'FORCE')
    
    #cc.add_parameter_to_local_subcase(subcase0_id, 'FORCE = ALL')
    cc.add_parameter_to_local_subcase(1, 'LOCAL_VAR = ALL')
    cc.add_parameter_to_local_subcase(1, 'LOCAL_VAR_B = ALL')
    #cc.add_parameter_to_global_subcase('GLOBAL_VAR = ALL')
    
    model_main.write_bdf('junk.bdf')
    sys.exit()

    # ------------------------------------------------------------------------
    # update the main bdf with the INCLUDE fiel
    
    contact_include = write_include_path(contact_bdf)
    model_main.rejects.append(contact_include)
    # ------------------------------------------------------------------------

    eid_start = self.get_max_eid() + 1 # elements + 1 -> starting ID
    eid = eid_start

    # ------------------------------------------------------------------------
    # apply springs based on distance

    f = open(contact_bdf, 'wb')
    eid_groups = []
    nodes_groups = []
    neids = []
    for contact_surface in contact_surfaces:
        (nids_group1, nids_group2, dof, stiffness, cid, glue, initial_gap, max_deflection_error) = contact_surface

        #model_left = BDF()
        #model_left.read_bdf(left_bdf)

        #model_right = BDF()
        #model_right.read_bdf(right_bdf)

        # change nodes from cid=0 to cid=N
        nodes_left  = update_nodes_cid(model_main, nids_group1, cid)
        nodes_right = update_nodes_cid(model_main, nids_group2, cid)

        # find 5 closest nodes; TODO: update this...
        spring_sets = {
        # left  # right
            1 : [2, 4, 5],
            2 : [10, 3],
        }

        c1 = c2 = dof
        neids_start = eid
        i = 0
        for g1, spring_nodes in sorted(spring_sets.iteritems()):
            f.write('$ contact set %i' % i)
            for g2 in spring_nodes:
                celas = ['CELAS2', eid, stiffness, g1, c1, g2, c2]
                f.write(print_card(celas))
            eid += 1
            i += 1
        neid = eid - neids_start
        
        eid_group = (spring_eids_min, spring_eids_max)
        eid_groups.append(eid_group)
        nodes_group = (nodes_left, nodes_right)
        nodes_groups.append(nodes_group)
        neids.append(neid)
    f.close()

    # ------------------------------------------------------------------------
    # we put this at the bottom b/c we need to update the node coordinate systems
    model_main.write_bdf(main_bdf2)
    del model_main
    
    # ------------------------------------------------------------------------
    # set the initial stiffnesses
    stiffnesses = zeros(sum(neids), 'float64')
    errors = ones(sum(neids), 'float64')
    ieid = 0
    for neid, contact_surface in zip(neids, contact_surfaces):
        left_bdf, right_bdf, dof, stiffness, cid, glue, initial_gap = contact_surface
        stiffnesses[ieid : ieid + neid] = stiffness
        ieid += neid
    return eid_groups, nodes_groups, stiffnesses, errors

def parse_op2(contact_bdf, op2_name, subcase_id, contact_surfaces, eid_groups, nodes_groups, stiffnesses, errors):
    """
    assumes static analysis
    """
    op2 = OP2()
    op2.read_op2(op2_name)
    spring_forces = op2.spring_forces[subcase_id]
    
    #force_name = 'forces.bdf'

    contact_bdf = 'spring.bdf'
    f = open(contact_bdf, 'wb')

    errorsi = 0
    for contact_surface, eid_group, nodes_group in zip(contact_surfaces, eid_groups, nodes_groups):
        (left_bdf, right_bdf, dof, stiffness, cid, glue, initial_gap, max_deflection_error) = contact_surface
        spring_eids_min, spring_eids_max = eid_groups
        nodes_left, nodes_right = nodes_groups

        nelements = spring_eids_max - spring_eids_min
        #forces = zeros((nelements), 'float64')
        c1 = c2 = cid

        i = 0
        for eid, force in spring_forces.forces.iteritems():
            #forces[i] = force
            g1 = spring_forces.g1[eid]
            g2 = spring_forces.g2[eid]

            
            k = stiffnesses[i]
            deflection = force/k
            if force <= 0:
                # assume a really small stiffness (no contact)
                #force = 0.0
                #if errors[i] == 1:  # if old error
                new_stiffness = 0.01
                stiffnesses[i] = new_stiffness
                new_flag = 0
            else:  # force > 0
                # the contact is correct
                #force_card = ['FORCE']
                #if deflection > 

                new_stiffness = stiffness
                #flag = int(force/abs(force)) # 0, 1
                new_flag = 1

            celas = ['CELAS2', eid, new_stiffness, g1, c1, g2, c2]
            f.write(print_card(celas))

            # check to see
            if new_flag != errors[i]:
                errorsi += 1
                errors[i] = new_flag
            i += 1
    f.close()
    ierrrors = where(errors == 1)[0]
    nerrors = len(ierrors)
    return nerrors

if __name__ == '__main__':
    model = BDF()
    model.read_bdf('plate.bdf')
    model.write_bdf('plate2.bdf')

    nids = [19, 20, 21, 22, 23, 24,]
    nids = model.nodes.keys()
    nids.sort()

    def centroid_at_7(element):
        c = element.Centroid()
        if c[1] == 7.0:  # update if the element centroid > 6
            #print c
            return True  # upper group
        return False     # lower group
    func = centroid_at_7
    nnodes = len(model.nodes)

    if 0:
        split_model(model, nids, func)
        model.write_bdf('plate_split.bdf')

        model2 = BDF()
        model2.read_bdf('plate_split.bdf')
        assert nnodes + 12 == len(model2.nodes), 'nnodes+12=%s nnodes2=%s' % (nnodes + 12, len(model2.nodes))
        #assert nnodes + 6 == len(model2.nodes), 'nnodes+6=%s nnodes2=%s' % (nnodes + 6, len(model2.nodes))
    
    group1_nodes = set([])
    group2_nodes = set([])
    for eid, element in model.elements.iteritems():
        nids = element.nodeIDs()
        if func(element):
            for nid in nids:
                group1_nodes.add(nid)
        else:
            for nid in nids:
                group2_nodes.add(nid)

    group1_nodes = list(group1_nodes)
    group2_nodes = list(group2_nodes)
    surface1 = {
        'group1_nodes' : group1_nodes,
        'group2_nodes' : group2_nodes,
        'dof'          : 1,     # dof in contact (1, 2, 3, 4, 5, 6) - nodes in group1/2
        'K'            : 1.e8,  # stiffness of contact interface
        'cid'          : 1,     # direction of contact
        'glue'         : False, # 
        'initial_gap'  : 0.0,
        'max_deflection_error' : max_deflection_error,
    }
    contact_surfaces = [surface1]
    main_bdf = 'plate_split.bdf'
    main_op2 = 'plate_split.op2'

    run(contact_surfaces, main_bdf, main_op2)