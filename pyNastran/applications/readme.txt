These applications are not maintained and are in various states of completion.
They are described in order of how current/stable they are:

 1. rainflow
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
