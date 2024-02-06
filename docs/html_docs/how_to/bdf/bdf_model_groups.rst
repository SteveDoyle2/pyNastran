Comments
--------
The standard nastran comment is supported.  Nastran does nothing with it.
```
$ CG
GRID        1001         742.959    270. 89.4568
```

Let's access that with code:
```python
bdf_model.nodes[1001].comment
```

Comments can be modified and written out.  Leave off the ``$``.

Model Groups
------------
Model groups intended as way to create tags for the GUI.  They can be used for other tasks though.

They're loaded into ``bdf_model.model_groups`` and left as a comment on the card.
As such, they're not written out because they'd lose their association with the data.

Some features:
 - An object may be in multiple groups at once
 - Element groups are recognized by ``--groups`` flag
 - node groups are recognized in the "Edit Geometry Properties" menu
 
Group types include:
 - nodes
 - elements
 - rigid_elements
 - masses
 - properties
 - materials
 - mpcs
 - spcs
 - loads

Correct:
```
$ group: name='RLongeron MainFuseStruct Gridpoint'; nodes=167:205
$ group: name='Skin MainFuseStruc'; elements=6001:15017
$ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123"; spcs=3
```

Incorrect due to a semicolon in the name:
```
$ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints; 123"; spcs=3
```

Unioned Groups (spcs becomes [1,2,3] internally):
```
$ group: name="spcs"; spcs=1
$ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, all"; spcs=1
SPC1           1  123456      13
$ group: name="spcs"; spcs=2
$ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123456"; spcs=2
SPC1           2  123456      13
$ group: name="spcs"; spcs=3
$ group: name="ULFuseCanardAtch MainFuseStruct Fixed point constraints, 123"; spcs=3
SPC1           3  123         13
```
