==========================
MSC/NX Nastran Differences
==========================

The following is a list of differences between MSC and NX Nastran.

 - OFFT field on the CBAR/CBEAM.  The OFFT field (default='GGG' for MSC)
   must be blank.  I'm pretty sure GGG is the only form NX supports.
 
 - In NX, the nodes that MPCs/RBEs must actually exist and have associated mass/stiffness.
 
 - In NX, massless nodes can mess results up, e.g., two springs attached at a massless grid.

