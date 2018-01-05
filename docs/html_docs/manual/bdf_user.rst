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

Attributes:
=================== ===========================================================
Name                 Description
=================== ===========================================================
``nid``             node id
``xyz``             Raw location <:math:`x_1, x_2, x_3`> in the BDF
``cp``              reference coordinate system
``cd``              analysis coordinate system
``ps``              nodal-based constraints
``seid``            superelement id
``cp_ref``          cross-referenced cp (or None)
``cd_ref``          cross-referenced cd (or None)
=================== ===========================================================

Attribute Methods:
=================== ===========================================================
Name                 Description
=================== ===========================================================
``xyz``             Raw location <:math:`x_1, x_2, x_3`> in the BDF
``Nid()``           gets nid
``Cp()``            gets cp_ref.cid or cp depending on cross-referencing
``Cd()``            gets cd_ref.cid or cd depending on cross-referencing
``Ps()``            gets ps
``SEid()``          superelement id
=================== ===========================================================


Methods:
==================== ===========================================================
Name                 Description
==================== ===========================================================
``get_position``     location in the global coordinate system
``get_position_wrt`` location in an arbitrary coordinate system
``set_position``     change the value of xyz and cp
``resolveGrids``     update all the positions of the grids to a common
                     coordinate system in an more efficient way
==================== ===========================================================


Using the GRID object::

 model = BDF()
 model.read_bdf(bdf_filename)
 node = model.Node(nid)
 node.get_position()                     % gets the position of the node in the global frame
 node.get_position_wrt(model, cid=0)     % same
 node.get_position_wrt(model, cid=1)     % gets the position of the node in a local frame

 node.set_position(model, array([1.,2.,3.]), cid=3) % change the location of the node

 model2 = BDF()
 model2.read_bdf(bdf_filename2)


CORDx
^^^^^^

A coordinate system may be defined by 3 points in any non-circular coordinate
system or by 3 vectors.

Once cross-referenced a node may use it's ``node.get_position_wrt(model, cid)`` method to
determine the position in any frame defined by a coordinate system.

