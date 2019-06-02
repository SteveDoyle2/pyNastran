if __name__ == '__main__':  # pragma: no cover
    data = """
/com,*********** Create Remote Point "Internal Remote Point 39" ***********
! -------- Remote Point Used by "Fixed - Line Body To EndCap 14054021-1 d" --------
*set,_npilot,803315
_npilot474=_npilot
et,332,170
type,332
real,332
mat,332
keyo,332,2,1              ! don't fix pilot node
keyo,332,4,0              ! MPC for all DOF's
tshape,pilo
en,501901,803315        ! create pilot node for rigid link
tshape
en,501902,803315,127827
/com,*********** Create Remote Point "Internal Remote Point 40" ***********
! -------- Remote Point Used by "Fixed - Line Body To EndCap 14054021-1 d" --------
*set,tid,334
*set,cid,333
et,cid,175
et,tid,170
keyo,tid,2,1               ! Don't fix the pilot node
keyo,tid,4,111111
keyo,cid,12,5              ! Bonded Contact
keyo,cid,4,0               ! Rigid CERIG style load
keyo,cid,2,2               ! MPC style contact
mat,333
real,333
type,333
en,501903,418114
en,501904,418115
en,501905,418116
en,501906,418117
en,501907,418118
en,501908,418119
en,501909,418120
en,501910,418121
en,501911,418122
en,501912,418123
en,501913,418124
en,501914,427511
en,501915,427512
en,501916,427518
en,501917,427524
en,501918,427528
en,501919,427533
en,501920,427539
en,501921,427544
en,501922,427551
en,501923,427562
en,501924,427569
*set,_npilot,803316
_npilot475=_npilot
type,tid
mat ,cid
real,cid
tshape,pilo
en,501925,_npilot
tshape



et,2,187
et,27,187   # element, group 27, element_type=187 -> tet10
et,30,188

etype   nastran_name
187     tet10
186     hexa20
188     beam

eblock,19,solid,,213
eblock,19,solid,,8
#----------------------------------------------------------------
et,_jid,184
et,tid,170
et,cid,174


keyo,tid,2,1               ! Don't fix the pilot node
keyo,tid,4,111111
keyo,cid,12,5              ! Bonded Contact
keyo,cid,4,2               ! Rigid CERIG style load
keyo,cid,2,2               ! MPC style contact
eblock,10,,,16

"""

from pyNastran.converters.dev.ansys.ansys import read_ansys
