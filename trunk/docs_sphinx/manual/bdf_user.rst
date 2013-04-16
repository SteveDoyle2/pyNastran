==========
User: BDF
==========


This document is intended to be a reference guide for users.


Nodes
------

GRID
^^^^^

In NASTRAN, nodes can be in the global or a local coordinate frame (rectangular,
cylindrical, spherical).  Each node may reference a different coordinate frame 
:math:`cp` for the reference coordinate frame (what the value is in the
:math:`x_1, x_2`, and :math:`x_3` fields) and for the analysis coordinate frame
:math:`aid`.
     
Node Specific Variable Names:

=================== ===========================================================
Name                 Description
=================== ===========================================================
``xyz``             Raw location <:math:`x_1, x_2, x_3`> in the BDF
``cp``              reference coordinate system
``cd``              analysis  coordinate system
``Position``        location in an arbitrary coordinate system
``UpdatePosition``  change the value of xyz and cp
``resolveGrids``    update all the positions of the grids to a common 
                    coordinate system in an more efficient way
=================== ===========================================================


Using the node object::

 bdf = BDF()
 bdf.readBDF(bdfName)
 node = bdf.Node(nid)
 node.Position()              % gets the position of the node in the global frame
 node.Position(cid=0)         % same
 node.Position(cid=1)         % gets the position of the node in a local frame
 
 node.UpdatePosition(bdf,array([1.,2.,3.]),cid=3) % change the location of the node
 
 bdf.resolveGrids(cid=0)      % change the xyz of all nodes to the same coordinate system as cid
 
 bdf2 = BDF()
 bdf2.readBDF(bdfNameAlt)
 bdf.unresolveGrids(bdf2)     % change the coordinte system back to the coordinate system in bdf2


CORDx
^^^^^^

A coordinate system may be defined by 3 points in any non-circular coordinate 
system or by 3 vectors.
     
Once cross-referenced a node may use it's ``node.Position(cid)`` method to 
determine the position in any frame defined by a coordinate system.

