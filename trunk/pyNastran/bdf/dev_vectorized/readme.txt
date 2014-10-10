This is an attempt at vectorization of the BDF class and is far from being complete.

Approach
==========
 - model.element (Element object) stores elements/properties/reference to model
   - Element
     - Element slice methods to get eids 1,2,3 (model.elements[3],
       model.elements[[1,2,3]] and model.elements[1:10:2]) never fail.
       They return None if values don't exist
     - model.get_elements will complain if results are not found
   - Mass
     - model.elements.get_mass([1,2,3]) does not fail
     - model.get_mass([1,2,3]) can fail

  - Easily vectorizable and/or high payoff cards (e.g. GRID, CTRIA3, CQUAD4, CTETRA4, CELAS1, PSHELL, PSOLID) are vectorized
    - When cards are grouped (e.g. CTRIA3, CQUAD4 are shells), a controlling class can be used to interface with the group of cards.
  
  - Complicated cards or low payoff cards are not (e.g. CORD2R, EIGRL, MAT1, MATS1) are not vectorized.
    The already implemented unvectorized BDF card will be used.

  - No cross referencing is allowed.  This may change for obscure cards.
