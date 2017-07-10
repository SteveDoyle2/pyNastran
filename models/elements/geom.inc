SPOINT,100,THRU,102
GRID     13             .5      .5      3.
$
GRID     5              1.      1.      1.
GRID     6               0.     1.      1.
GRID     7               0.      0.     1.
GRID     8              1.       0.     1.
$
GRID     12              0.      0.     2.
GRID     10              0.     1.      2.
GRID     9              1.      1.      2.
GRID     11             1.       0.     2.
$
GRID     1               0.      0.      0.
GRID     4               0.     1.       0.
GRID     3              1.      1.       0.
GRID     2              1.       0.      0.
$
GRID     14              0.      0.     -1.
GRID     15              0.     1.      -1.
GRID     16             1.      1.      -1.
GRID     17             1.       0.     -1.
$
GRID     18              0.      0.     -2.
GRID     19              0.     1.      -2.
GRID     20             1.      1.      -2.
GRID     21             1.       0.     -2.
$
GRID     22              0.      0.     -3.
GRID     23              0.     1.      -3.
GRID     24             1.      1.      -3.
GRID     25             1.       0.     -3.

GRID     26              0.      0.     -4.
GRID     27              0.     1.      -4.
GRID     28             1.      1.      -4.
GRID     29             1.       0.     -4.

CHEXA    1       2       2       3       4       1       8       5
         6       7
CPENTA   2       2       6       8       5       10      11      9
CPENTA   3       2       6       7       8       10      12      11
CTETRA   4       2       10      11      9       13
CTETRA   5       2       10      12      11      13
$ pshells
CQUAD4   6       4       4       1       14      15
CQUAD4   7       4       3       2       17      16
CTRIA3   8       4       4       3       16
CTRIA3   9       4       16      15      4
CTRIA3   10      4       1       2       17
CTRIA3   11      4       17      14      1
$
CBEAM    12      5       22      26      0.     1.       0.
CBAR     13      1       23      27      0.     1.       0.

$ rods
CROD     14      3       24      28
CROD     15      3       25      29
PROD     3       1      .1      .0088021 1.

$ pcomp
CQUAD4   16      6      14      15       19      18
CQUAD4   17      7      17      16       20      21
CTRIA3   18      6      18      14       17
CTRIA3   19      6      18      17       21
CTRIA3   20      7      19      15       16
CTRIA3   21      7      19      16       20
PCOMP,6,,,,,,,
,1,0.1
,1,0.2
,1,0.3
,1,0.4
PCOMP,7,,,,,,,
,1,0.1,,,1,0.2
,1,0.3,,,1,0.4
,1,0.5

$ shear
CSHEAR	22	8	19	20	24	23
PSHEAR	8	1	0.1

$ pcompg
CQUAD4	23	9	18	21	25	22
PCOMPG,9,,,,,,,
,1,1,0.1
,2,1,0.2
,4,1,0.3
,3,1,0.4

$ drilling
CQUADR	24	9	18	19	23	22
CTRIAR	25	9	21	25	20

$ conrod
CONROD	26	22	30	1	1.0	2.0
GRID     30             -1.      0.     -3.


$ ctube
GRID     31              2.      0.     -3.
CTUBE	28	10	25	31
PTUBE	10	1	0.1

$ celas
CELAS1	30	11	25	2	31	2
CELAS1	31	11	25	3	31	3
CELAS1	32	11	25	5	31	5
CELAS1	33	11	25	6	31	6
PELAS	11	100.
CELAS2	34	100.	22	6	30	6

CDAMP1	40	12	25	2	31	2
CDAMP1	41	12	25	3	31	3
CDAMP1	42	12	25	5	31	5
CDAMP1	43	12	25	6	31	6
PDAMP	12	2.
CDAMP2	44	12.	22	6	30	6

CDAMP3	45	13	101	102
CELAS3	46	14	101	102
CDAMP4	47	15.	102
CELAS4	48	16.	102
CELAS2	49	17.	25	1	101	0
PDAMP	13	13.
PELAS	14	14.

$CDAMP1	41	12	25	3	31	3
$CDAMP1	42	12	25	5	31	5
$CDAMP1	43	12	25	6	31	6
$PDAMP	12	2.
$CDAMP2	44	12.	22	6	30	6


$ visc
GRID     32              -1.	0.	-2.
GRID     33              2.	0.	-2.
CVISC	50	50	18	32
CVISC	51	51	21	33
PVISC	50	1.0	2.0		51	3.0	4.0
CONM2	50	32		1.0
CONM2	51	33		1.0

$ cquad8/ctria6

$  4--- x----60
$  |         |
$  |         62
$  |         |
$  15---63---61

$	eid	pid	n1	n2	n3	n4	n5	n6
$	n7	n8
CQUAD8	60	4	4	15	61	60		63
	62	65

$	eid	pid	n1	n2	n3	n4	n5	n6
CTRIA6	61	4	6	4	60		65	64
$GRID     4               0.     1.       0.
$GRID     15              0.     1.      -1.
GRID	60		-1.	1.	0.
GRID	61		-1.	1.	-1.
GRID	62		-1.	1.	-0.5
$ 15-62
GRID	63		-0.5	1.	-1.
GRID	64		-0.5	1.	0.5
GRID	65		-0.5	1.	-0.

$ bend
$CBEND	33	12	21	31	0.	0.	1.	1
$PBEND	12	1	1.0	1.0	0.0	0.0

PBEAM          5       1     .25.0052083.0052083        .0087949
            -.25    -.25     .25    -.25     .25     .25    -.25     .25
             YES      1.     .25 .005208 .005208      0..0087949      0.
        .0087949      0.    -.25    -.25     .25    -.25     .25     .25
         .850116.8501160      0.      0.      0.      0.      0.      0.
              0.      0.      0.      0.      0.      0.      0.      0.
PBAR     1       1      0.25    .0052083.0052083.0088021
PSOLID   2       1       0
PSHELL   4       1      .25      1               1
