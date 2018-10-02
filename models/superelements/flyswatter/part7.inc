$
$ part7.dat
$
CQUAD4	1	1	1	2	4	3
CQUAD4	2	1	3	4	6	5
CQUAD4	3	1	5	6	8	7
CQUAD4	4	1	7	8	14	13
$
$		cp	a	b	c	cd	ps	seid
GRDSET							6
$
$	nid	cp	x	y	z	cd	ps
GRID	1		-.4	0.	0.		123456
GRID	2		.4	0.	0.		123456
GRID	3		-.4	0.9	0.
GRID	4		.4	0.9	0.
GRID	5		-.4	1.8	0.
GRID	6		.4	1.8	0.
GRID	7		-.4	2.7	0.
GRID	8		.4	2.7	0.
$
$ boundary grids
$
GRID	13		-.4	3.6	0.
GRID	14		.4	3.6	0.
$
MAT1,1,30.+6,,.3,.283
PARAM,WTMASS,.00259
PARAM,AUTOSPC,YES
PSHELL,1,1,.05,1,,1