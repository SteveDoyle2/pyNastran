
"""
animation
Case 1:  "Transient Result"

Results Sidebar
---------------
Subcase 1
 - Time 1
    - Disp 1
    - Stress 1
 - Time 2
    - Disp 2
    - Stress 2


Cases Added
-----------
Stress 1
Stress 2
...
Stress N


Add Cases From:  Subcase 1
Add Result:      Disp
Mode:            Transient
Method:          Actual/DispScale=40


Add Cases From:  Subcase 1
Add Result:      Stress
Mode:            Transient
Method:          Actual
Case 2:
=======

Results Sidebar
---------------
Subcase 1
 - Time 1
    - Eigenvectors 1
Subcase 2
 - Frequency 1
     - Disp 1

Add Cases From:  Subcase 1
Add Result:      Stress 1
Mode:            Eigenvector
Method:          Cosine Ramp/Linear
"""

