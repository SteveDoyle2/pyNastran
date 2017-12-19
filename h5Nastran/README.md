h5Nastran is technically not a part of pyNastran (yet), but uses pyNastran as a dependency.  It is in the early stages of development.  Contributions are welcome!

h5Nastran is a python package that deals with the MSC.Nastran h5 file format (2018 version).  It uses pyNastran to convert bdf files to the h5 format, and has a standalone punch and f06 results reader to convert to the h5 format.  Results tables are easily searchable.

Example:
```python
from h5Nastran import H5Nastran

db = H5Nastran('./models/model_001.h5', 'w')
db.load_bdf('./models/model_001.bdf')
db.load_punch('./models/model_001.pch')

print(db.input.node.grid.identity)  # or db.input.node.grid.grid

domain_ids = [1, 2]
elements = [400002, 400111, 400198]
forces = db.result.elemental.element_force.quad4.search(domain_ids, elements)

print(forces)

# pynastran bdf
bdf = db.bdf

# currently, to modify the bdf and rewrite to h5,
# you'd need to modify the pynastran bdf, write it to a new file and create a new h5 database

# currently, the entire bdf is written to h5 as written by pynastran
# it can be loaded by doing db.load_bdf() without a filename
# the goal is to recreate the bdf file from only the h5 data

db.close()
```


Not all bdf features are currently supported.  Bulk data cards are pretty easy to add... either make a pull request or request one to be added.  Punch tables are also easy to add - either make a pull request or request one to be added (but please provide an example of the punch table format).

bdf's, punch files, and f06 files are appreciated so I can add more support and do more testing.  The punch file reader is currently much better supported than the f06 file reader.

There is no intention to support SORT2 results tables unless someone can make a very good case for supporting them.

Roadmap:
1.  add more bdf cards and other bdf information support (control decks, etc)
2.  add more resuls tables (punch files first priorty, f06 second priority)
3.  op2 to h5?  same idea would be used as for f06 reader
4.  add some more results search features
5.  create standalone GUI to query results (no VTK or plotting), just data
6.  once GUI is done, it can be incorporated into pyNastran GUI if desired
