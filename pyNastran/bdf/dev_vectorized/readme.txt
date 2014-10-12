This is an attempt at vectorization of the BDF class and is far from being complete.

Approach
==========
 - clean break from old BDF class
   - Why?
     - the data structure is fundamentally different
   - Note: the CaseControlDeck and ExecutiveControlDeck will be unchanged
 
 - model.element (Element object) stores elements/properties/reference to model
   - Element
     - Element slice methods to get eids 1,2,3 (model.elements[3],
       model.elements[[1,2,3]] and model.elements[1:10:2])
       For missing elements/properties, it returns:
         - Option A:  None if values don't exist
         - Option B:  No entry
         - Option C:  Crash
         - no decision has been made on this, but leaning towards Option B
         - Why A?
           - you don't have a valid BDF
             - nobody ever has a valid model...
           - don't waste developer time adding checks for getting the right
             number of properties that were requested (it'll probably crash on None)
           - no need to return eids when mass is called
         - Why B?
           - you should have a valid BDF / code
           - it should be a lot faster
           - should we return eids when mass is called?
         - Why C?
           - best option for developer time.  
           - a crash flag might be good
         - Note
           - It's been specifically requested that code doesn't crash on
             missing references and that None may be a valid option
           - You can look at geometry without looking at properties/materials
           - However, for something like Mass, maybe this doesn't make sense
           
     - model.get_elements will complain if exact list of elements are not found?
       - could remove this method if we figure A/B/C out properly
   - Mass
     - model.elements.get_mass([1,2,3]) does not fail???
     - model.get_mass([1,2,3]) can fail???
   - Fail Criteria
     - regardless of A/B/C; stills fails given very bad input
       - strings when values should be integers
       - None (when not explicitly allowed)
     - if total mass is desired:
       - should it fail for None values?
       - should it set those values to 0.0?
       - or just not sum the mass?
   - Why the differing fail criteria?
     - it requires extra checks; speed
     - if you as Patran/Femap to show elements 1-10 and only element 1 exists,
       it shows you element 1 and (hopefully :) doesn't crash

  - Easily vectorizable and/or high payoff cards
    - These include:
       - GRID, SPOINT (nodes)
       - CTRIA3, CQUAD4, CTRIA6, CQUAD8 (shells)
       - CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20 (solids)
       - CELAS1, CELAS2, CELAS3, CELAS4 (springs)
       - CONROD, CROD (rods)
       - PELAS
       - PSHELL, PCOMP, PCOMPG (shell properties)
       - PSOLID, PLSOLID (solid properties)
    - When cards are grouped (e.g. CTRIA3, CQUAD4 are shells),
      a controlling class is used to interface with the group of cards.
      - Why?
        - interfacing is a pain for vectors
        - different people want to interface differently
          - solver
          - pyNastran object style
  
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
    - Note:
       - methods on the original cards (if used) will be updated to work
         without cross referencing
       - APIs will be updated to be as similiar as possible when a separate card is used

  - No cross referencing is allowed.  This may be OK for obscure cards.
    - Why?
      - cross referencing prevents some options unless you re-cross-reference
        - renumbering
      - slow
      - it prevents vectorization
