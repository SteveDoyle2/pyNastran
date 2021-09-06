This is an attempt at vectorization of the BDF class and is far from being complete.

Approach
==========
 o clean break from old BDF class
   o Why?
     o the data structure is fundamentally different
   o Note: the CaseControlDeck and ExecutiveControlDeck will be unchanged
 
 o model.element (Element object) stores elements/properties/reference to model
   o Element
     o Goals...
       o Get elements by element id
         o model.elements[3]
         o model.elements[[1,2,3]]
         o model.elements[1:10:2])
       o Get elements by index
         o ???
       o semi sanely handle missing data
        
       o Handling missing data:
         o Possible Options:
           o Option A:  NaN if values don't exist
           o Option C:  Crash on Failure
           o Option D:  Crash Flag (A or C)
           o Option E:  Validate Inputs when model is loaded (as a flag)
         o Invalid Options:
           o Option B:  No entry
             o removed as a potential Option in D

         o a final decision has been made on this, but
           o old pyNastran BDF() uses Options C/E
           o Jeff has requested not C, forced E
             o possibly A/B/D/E
           o leaning towards Option E (with a GEOMCHECK/input flag)
             and D
           o B is NOT going to happen.  It has too many downsides.

         o Why A? (NaN value)
           o NaN could possibly be None
           o you don't have a valid BDF
             - nobody ever has a valid model...
           - wastes developer time adding checks for getting the right
             number of properties that were requested (it'll probably crash on NaN)
           + no need to return eids when mass is called
         o Why B? (No entry)
           o you should have a valid BDF / code
           + it's faster because:
             + elements/properties are eliminated and don't need to be handled later
           - should we return eids when mass is called?
           - errors pass silently
         o Why C? (Crash on Failure)
           o doesn't fail unless a bug is specifically encountered
           + allows for incomplete models
           - requires checks in more places,
           + which probably should exist anyways
           - checks are inevitably incomplete and poor error messages will result
         o Why D? (Validate Inputs + Crash Flag)
           + compromise between A and C
           o still need to decide on A vs. C
         o Why E? (Validate Inputs)
           + best option for avoiding bugs
           - could be tedious
           - Jeff has specifically requested against this option
           o a crash flag (e.g. GEOMCHECK) could be used (option D)

         o Results
           o Approved: Option E as a flag
           o Eliminated: Option B

           o Jeff has requested not C
             o possibly A
           o Conclusion
             o Option A/C as a flag is probably the way to go
             o since you gotta pick one to start with and Option A
               promotes bugs, Option C will be implemented first

           o request for elements that don't exist
             o use case #1:  user requests mass from elements 1:10, but 5 doesn't exist
                 o Favorite:   Option C
                 o Compromise: Option A
 
           o request for elements that don't exist
             o use case #2:  user requests mass from elements 1:10, but
             		     Node 1 / Material 1 on element 1 doesn't exist
                 o Favorite:   Option C
                 o Compromise: Option A

           o non-request for nodes/materials that don't exist
             o use case #3:  user loads model, but missing some data
                 o Favorite:   Option A/C

             o use case #4:  user wants to write elements 1:10, but 5 doesn't exist
                 o Favorite:   Option A
                 o Compromise: Option C

         o What about?
           o request for elements that don't exist
             o use case #1:  user requests mass from elements 1:10, but 5 doesn't exist
                 o Favorite:   Option C/E
                 o Compromise: Option A
                 o Eliminated: Option B

                 o Why #1: intentionally (e.g. element 1:#)
                    solution A: NaN value
                      + they need to deal with the value explicitly, but got what
                        they wanted (good enough)
                    solution B: drop the element -> No entry
                      + they got what they wanted (ideal)
                    solution C: crash on failure
                      - they have to fix it (annoying)
                    Favorite:   Option B
                    Compromise: Option A
                    Dislike:    Option C

                 o Why #2: they have a bug
                    solution A: NaN value
                      + they can see the NaN, but could get confused with Why #1
                        (good enough)
                    solution B: drop the element -> No entry
                      - errors pass silently (bad!)
                      - off by 1 errors WILL happen, which forces all sorts of checks
                         which makes code complicated, hard to maintain, and buggy
                         (bad!!!)
                    solution C/E: crash on failure / validate input
                      + the error is caught with an error message (ideal)
                    Favorite:   Option C/E
                    Compromise: Option A
                    Eliminated: Option B

           o request for elements that don't exist
             o use case #2:  user requests mass from elements 1:10, but
             		     Node 1 / Material 1 on element 1 doesn't exist
                 o Favorite:   Option C/E
                 o Compromise: Option A
                 o Eliminated: Option B
                 o Why #1: bug
                    solution A: NaN value
                      + they get a NaN (good enough)
                    solution B: drop the element -> No entry
                      - errors pass silently (bad!!!)
                    solution C: crash on failure
                      + the error upon requesting mass and is caught
                        with an error message (ideal)
                    solution E: validate inputs
                      + the user never gets past the loading step
                        and gets an error message (ideal)
                    Favorite:   Option C/E
                    Compromise: Option A
                    Eliminated: Option B

     
           o non-request for nodes/materials that don't exist
             o use case #3:  user loads model, but missing some data
                 o Favorite:   Option A/C
                 o Compromise: Option E

                 o Why #1: Volume of solids is desired, so I don't need density
                           on the undefined MAT1 card
                    solution A/C: NaN value / crash on failure
                      + they got what they wanted (ideal)
                    solution E: validate inputs
                      - the user never gets past the loading step
                        and gets an error message (bad)
                    Favorite:   Option A/C
                    Compromise: Option E

                 o Why #2: Volume of solids is desired, so I don't need the location
                           of Node 1 used on the undefined CELAS1 card
                    solution A/C: NaN value / crash on failure
                      + they got what they wanted (ideal)
                    solution E: validate inputs
                      - the user never gets past the loading step
                        and gets an error message (bad)
                    Favorite:   Option A/C
                    Compromise: Option E

             o use case #4:  user wants to write elements 1:10, but 5 doesn't exist
                 o Favorite:   Option A
                 o Compromise: Option C
                 o Eliminated: Option B

                 o Why #1: intentionally (e.g. element 1:#)
                    solution A: NaN value
                      + print a flag such as "CQUAD4 5 doesnt exist" (good enough)
                    solution B: drop the element -> No entry
                      + they got what they wanted (ideal)
                    solution C: crash on failure
                      - they have to fix it (annoying)
                    Favorite:   Option B
                    Compromise: Option A
                    Dislike:    Option C

                 o Why #2: they have a bug
                    solution A: NaN value
                      + they can see the warning (good enough)
                    solution B: drop the element -> No entry
                      - errors pass silently (bad!)
                    solution C/E: crash on failure / validate input
                      + the error is caught with an error message (ideal)
                    Favorite:   Option C
                    Compromise: Option A
                    Eliminated: Option B

           
     o model.get_elements will complain if exact list of elements are not found?
       o could remove this method if we figure A/B/C/D/E out properly
   o Mass
     o model.elements.get_mass([1,2,3]) does not fail???
     o model.get_mass([1,2,3]) can fail???
   o Fail Criteria
     o regardless of A/B/C/D/E; stills fails given very bad input
       o strings when values should be integers
       o None (when not explicitly allowed)
     o if total mass is desired:
       o should it fail for None values?
       o should it set those values to 0.0?
         - sounds like Option B (errors pass silently)
       o or just not sum the mass?
   o Why the differing fail criteria?
     o it requires extra checks; speed
     o if you as Patran/Femap to show elements 1-10 and only element 1 exists,
       it shows you element 1 and (hopefully :) doesn't crash
   o Interface
     o 2D arrays will iterate over the rows to indicate different elements/properties
     o column 2/3/etc. will refer to column 1
       o yes [eid1, mass1]
             [eid2, mass2]
             [eid3, mass3]
             [eid4, mass4]
       o no  [eid1,   eid2,  eid3,  eid4]
             [mass1, mass2, mass3, mass4]
     o Why?
       o I'll forget which is which and a standard is a good thing

  o Easily vectorizable and/or high payoff cards
    o These include:
       o GRID, SPOINT (nodes)
       o CTRIA3, CQUAD4, CTRIA6, CQUAD8 (shells)
       o CTETRA4, CPENTA6, CHEXA8, CTETRA10, CPENTA15, CHEXA20 (solids)
       o CELAS1, CELAS2, CELAS3, CELAS4 (springs)
       o CONROD, CROD (rods)
       o PELAS
       o PSHELL, PCOMP, PCOMPG (shell properties)
       o PSOLID, PLSOLID (solid properties)
    o When cards are grouped (e.g. CTRIA3, CQUAD4 are shells),
      a controlling class is used to interface with the group of cards.
      o Why?
        + interfacing is a pain for vectors
        + different people want to interface differently
          o solver
          o pyNastran object style
        + it doesn't even matter as users can choose to not use these classes
  
  o Complicated cards or low payoff cards will not be vectorized.  The already
    implemented unvectorized BDF card will be used (ideally).
    o These include:
      o CORD2R, MAT1, MATS1
      o CBEAM
      o PCOMP, PCOMPG, PBEAML, PBARL
      o AEROS, CAERO1, EIGRL
    o Why?
       o Some of these cards are dynamic in length (PCOMP, CBEAM)
       o Some are interlinked and require special code (MAT1, MATS1)
       o Others are just a waste of time because you don't use that many
         (e.g. CORD2R, AEROS, CAERO1, EIGRL)
    o Note:
       o methods on the original cards (if used) will be updated to work
         without cross referencing
       o APIs will be updated to be as similar as possible when a separate
         card is used

  o No cross referencing is allowed.  This may be OK for obscure cards.
    o Why?
      o cross referencing prevents some options unless you re-cross-reference
        - renumbering
        - cards added after cross referencing
      - cross referencing cards is slow
      - accessing data for cross referenced cards is slow
      - writing out cross referenced cards is slow
      - it prevents vectorization



Vectorized Cards (done)
=======================
GRID

# mass
CONM1, CONM2

# Elements - 0D
PELAS, CELAS1, CELAS2, CELAS3, CELAS4

# Elements - 1D
PROD, CROD, CONROD
PBAR, CBAR
CBEAM
PBUSH, CBUSH

# Elements - 2D
PSHEAR, CSHEAR
PSHELL, CQUAD4, CTRIA3

# Elements - 3D
PSOLID, CTETRA4,  CPENTA6,  CHEXA8,
        CTETRA10, CPENTA15, CHEXA20

# Loads - 0D
FORCE, MOMENT, FORCE1, MOMENT1, GRAV

# Loads - 1D
PLOAD1, PLOAD2, PLOADX1, RFORCE
# Loads - 2D/3D

Vectorized Cards (not done)
===========================
# mass

# Elements - 1D
PDAMP, CDAMP1, CDAMP2, CDAMP3, CDAMP4, CDAMP5

# Elements - 2D
CTRIA6, CTRIAX6, CTRIA, CQUAD, CQUAD8, CQUADX

# Loads - 1D
# Loads - 2D/3D

Unvectorizable Cards (done)
===========================
# Elements - 1D
PBEAM, PBEAML, PBCOMP

# Elements - 2D
PCOMP, PCOMPG

Unvectorizable Cards (not done)
===============================
TABLEx

Vectorizable, but may not be worth vectorizing
==============================================
# aero
CAERO1, CAERO2, CAERO3, CAERO4, CAERO5
PAERO1, PAERO2, PAERO3, PAERO4, PAERO5

Not Grouped
===========
# mass
CMASS1, CMASS2, CMASS3, CMASS4, CMASS5

# aero
AERO, AEROS
SPLINE1, SPLINE2, SPLINE3, SPLINE4, SPLINE5
AEFACT, AESURF, AESURFS, AELIST, AEPARM, AESTAT, AELINK
FLFACT, GUST, MKAERO1, MKAERO2
TRIM, CSSCHD

 
# Elements - 1D
PDAMPT

# Elements - 2D
# Elements - 3D
# Loads - 1D
# Loads - 2D/3D

