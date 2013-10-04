from __future__ import division
import os

from numpy import zeros, ones, array, where

from pyNastran.bdf.fieldWriter import print_card
from pyNastran.bdf.writePath import write_include


max_deflection_error = 0.01
def run():

    group1_nodes = [1, 2, 3]
    group2_nodes = [4, 5, 6]
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
    contact_bdf = 'spring.bdf'
    main_bdf = 'main.bdf'
    op2_name = 'main.op2'
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

    subcases = model.case_control_deck.get_subcases()
    assert len(subcases) == 1, len(subcases)
    subcase0_id = subcases.keys()[0]
    subcase = subcases[subcase0_id]
    subcase['FORCE'] = ['FORCE', 'ALL']

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
        ]

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