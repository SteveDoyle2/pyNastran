These applications are not maintained and are in various states of completion.
They are described in order of how current/stable they are:

 1. aero_panel_buckling (created 11/2016)
     - Very complete & validated; requires a few parameters to setup a new model, but very general
     
     - Takes a BDF/OP2 of a static aero (144) solution and with minimal input creates nskin_panels based on:
      - geometry boundaries
      - symmetry plane (optional)
      - property boundary (optional)
     - Uses a very robust 2D geometry splitting algorithm, determines skin panels.
     - Extracts the deflection for each skin panel from the OP2 and generates N buckling BDFs
       with applied deflections and dummy loads
     - Determines the critical buckling eigenvalue on a per patch basis
    - needs a public example problem
    - created during v0.8 development
    - ignores solid & line elements
    
 2. rainflow
     - Very complete and validated
     - rainflow_from_csv
       - performs rainflow counting per ASTM E1049 5.4.5.2 on a generic stress/load spectrum
       - creates figures
     - no numpy array option, but it should be simple to support that
    - created during v0.7 development
    - enhanced during v0.9; added tests

 2. hyper
    - don't really remember how to run it, but it does work
    - implemented Newtonian impact -> Nastran PLOAD cards
    - needs a public example problem
    - created around v0.7-v0.8

 3. cart3d_nastran_fsi
    - This was the project that stated pyNastran
    - Used to be fully implemented and validated, but needs some tests
    - Likely broken code to perform aero mapping from Cart3d to Nastran
    - Runs Nastran
    - Maps structural deflections to Cart3d
    - Iterates to convergence
    - needs a public example problem
    - created around v0.1 (updated since)
    - updated load mapping part to be far more general
