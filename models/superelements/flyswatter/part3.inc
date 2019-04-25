$
$ part3.dat
$
CQUAD4	14	1	19	20	30	29
CQUAD4	15	1	29	30	36	35
$
$		cp	a	b	c	cd	ps	seid
GRDSET							6
$
$ boundary grids
$
GRID	19		-3.6	4.4	0.
GRID	20		-2.8	4.4	0.
GRID	35		-3.6	6.	0.
GRID	36		-2.8	6.	0.
$
GRID	29		-3.6	5.2	0.
GRID	30		-2.8	5.2	0.
$
MAT1,1,30.+6,,.3,.283
PARAM,WTMASS,.00259
PARAM,AUTOSPC,YES
PSHELL,1,1,.05,1,,1
$