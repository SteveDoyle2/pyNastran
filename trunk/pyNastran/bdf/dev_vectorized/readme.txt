This is an attempt at vectorization of the BDF class and is far from being complete.

Approach
==========
 - clean break from old BDF class
   - Why?
     - the class is fundamentally different
 
 - model.element (Element object) stores elements/properties/reference to model
   - Element
     - Element slice methods to get eids 1,2,3 (model.elements[3],
       model.elements[[1,2,3]] and model.elements[1:10:2]) rarely fails
       They return None if values don't exist
     - model.get_elements will complain if results are not found
   - Mass
     - model.elements.get_mass([1,2,3]) does not fail
     - model.get_mass([1,2,3]) can fail
   - Fail Criteria
     - should model fail given very bad input?
       - strings where values should be integers
       - None
     - if total mass is desired:
       - should it fail for None values
       - should it set those values to 0.0
       - or just not sum the mass?
   - Why the differing fail criteria?
     - it requires extra checks; speed
     - if you as Patran to show elements 1-10 and only element 1 exists, it shows you element 1

  - Easily vectorizable and/or high payoff cards
    - These include:
       - GRID, SPOINT
       - CTRIA3, CQUAD4
       - CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20
       - CELAS1, CELAS2, CELAS3, CELAS4, CONROD, CROD
       - PELAS, PSHELL, PSOLID
    - When cards are grouped (e.g. CTRIA3, CQUAD4 are shells),
      a controlling class can be used to interface with the group of cards.
      - Why?
        - interfacing is a pain for vectors
        - different people want to interface differently (e.g. solver vs. pyNastran object style)
  
  - Complicated cards or low payoff cards will not be vectorized.  The already
    implemented unvectorized BDF card will be used (ideally).
    - These include:
      - CORD2R, MAT1, MATS1
      - CBEAM
      - PCOMP, PCOMPG, PBEAML, PBARL
      - AEROS, CAERO1, EIGRL
    - Why?
       - Some of these cards are dynamic in length (PCOMP, CBEAM)
       - Some are interlinked and require special code (MAT1, MATS1)
       - Others are just a waste of time because you don't use that many
         (e.g. CORD2R, AEROS, CAERO1, EIGRL)

  - No cross referencing is allowed.  This may change for obscure cards.
    - Why?
      - cross referencing prevents some options unless you re-cross-reference
      - it prevents vectorization
