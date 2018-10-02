CQUAD4	5	1	13	14	24	23
$
$		cp	a	b	c	cd	ps	seid
GRDSET							6
GRID	13		-.4	3.6	0.
GRID	14		.4	3.6	0.
GRID	23		-.4	4.4	0.
GRID	24		.4	4.4	0.
$
MAT1,1,30.+6,,.3,.283
PSHELL,1,1,.05,1,,1
$
param,wtmass,.002588
