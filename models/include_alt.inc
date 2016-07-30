FORCE          1       3       0      1.    100.    200.    300.
SPC1           1  123456       1
SPC1           1  123456       2
PCOMP          1      0.      0.                                        +       
+              1      .1      0.     YES      10      .5      0.     YES+       
+              1      .1      0.     YES
MAT8           1    1.+7    1.+7      .33840000.3840000.3840000.      .1+       
+          2.4-5   2.4-5      0.                                        +       
+                               
MAT8          10      1.      1.      0.      1.  85000.   2300.      .2+       
+           5.-5    5.-5      0.                                        +       
+                               
$ PCOMP test
CQUAD4         1       1       1       2       3       4      0.
CTRIA3         2       1       1       2       3

$ PSHELL, NSM = 0.0
PSHELL         2       1    .125       1               1              0.
CQUAD4         3       2       1       2       3       4      0.
CTRIA3         4       2       1       2       3

$ PSHELL, NSM = 1.0
PSHELL         3       1    .125       1               1              1.
CQUAD4         5       3       1       2       3       4      0.
CTRIA3         6       3       1       2       3

$ Solid Nodes
GRID          11      0      0.      0.      0.       0
GRID          12      0      1.      0.      0.       0
GRID          13      0      1.      1.      0.       0
GRID          14      0      0.      1.      0.       0

GRID          15      0      0.      0.      2.       0
GRID          16      0      1.      0.      2.       0
GRID          17      0      1.      1.      2.       0
GRID          18      0      0.      1.      2.       0

$ Solids
PSOLID   4       1       0
CHEXA          7       4      11      12      13      14      15      16
              17      18
CTETRA         8       4      11      12      13      15


$ Solid Nodes
GRID          21      0      0.      0.      0.       0
GRID          22      0      1.      0.      0.       0
GRID          23      0      1.      1.      0.       0
GRID          24      0      0.      0.      2.       0
GRID          25      0      1.      0.      2.       0
GRID          26      0      1.      1.      2.       0
CPENTA         9       4      21      22      23      24      25      26
INCLUDE 'include2.inc'