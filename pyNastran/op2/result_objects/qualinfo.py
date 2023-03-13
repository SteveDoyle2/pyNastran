from typing import Any


class QualInfo:
    def __init__(self,
                 AUXMID=0,
                 AFPMID=0,
                 DESITER=0,
                 HIGHQUAL=0,
                 DESINC=0,
                 MASSID=0,
                 ARBMID=0,
                 TRIMID=0,
                 MODULE=0,
                 FLXBDYID=0,
                 ISOLAPP=0,
                 #ISOLAPP=1,
                 SEID=0,
                 PEID=0,

                 DISCRETE=False,
                 PRESEQP=False,
                 DELTA=False,
                 BNDSHP=False,
                 ADJOINT=False,
                 APRCH=' ',
                 PARTNAME=' ',
                 DFPHASE=' ',
                 QCPLD=' ',
                 DBLOCFRM=' ',
                 P2G=' ',
                 K2GG=' ',
                 M2GG=' ',
                 CASEF06=' ',

                 NL99=0,
                 MTEMP=0,
                 SUBCID=0,
                 #OSUBID=1,
                 OSUBID=0,
                 STEPID=0,
                 RGYRO=0,
                 SSTEPID=0,
                 #----------------------------------
                 PRELOAD=0,
                 MPC=0,
                 SPC=0,
                 SUPORT=0,
                 BMETH=0,
                 NLOOP=0,
                 STATSUB=0,
                 ):
        # ints
        self.AUXMID = AUXMID
        self.AFPMID = AFPMID
        self.DESITER = DESITER
        self.HIGHQUAL = HIGHQUAL
        self.DESINC = DESINC
        self.MASSID = MASSID
        self.ARBMID = ARBMID
        self.TRIMID = TRIMID
        self.MODULE = MODULE

        # booleans
        self.DISCRETE = DISCRETE
        self.PRESEQP = PRESEQP
        self.DELTA = DELTA
        self.BNDSHP = BNDSHP
        self.ADJOINT = ADJOINT

        # strings
        self.PARTNAME = PARTNAME
        self.DFPHASE = DFPHASE
        self.APRCH = APRCH
        self.QCPLD = QCPLD
        self.DBLOCFRM = DBLOCFRM
        self.P2G = P2G
        self.K2GG = K2GG
        self.M2GG = M2GG
        self.CASEF06 = CASEF06

        # other
        self.FLXBDYID = FLXBDYID
        self.ISOLAPP = ISOLAPP
        self.SEID = SEID
        self.PEID = PEID
        self.NL99 = NL99
        self.MTEMP = MTEMP
        self.SUBCID = SUBCID
        self.OSUBID = OSUBID
        self.STEPID = STEPID
        self.RGYRO = RGYRO
        self.SSTEPID = SSTEPID

    @classmethod
    def from_str(self, db_key: int, qual_str: str, log):
        assert qual_str[0] == '(', f'qual_str={qual_str!r}'
        assert qual_str[-1] == ')', f'qual_str={qual_str!r}'
        qual_str2 = qual_str[1:-1]
        slines = qual_str2.split(';')

        def _parse(slines: list[str], debug: bool=False) -> dict[str, Any]:
            data_dict = {}
            for i, sline in enumerate(slines):
                key, value_str = sline.split('=', 1)
                if value_str == 'FALSE':
                    value = False
                elif value_str == "TRUE":
                    value = True
                elif value_str == "' '":
                    value = ''
                else:
                    value_str = value_str.strip()
                    if value_str.isdigit():
                        value = int(value_str)
                    else:
                        value = value_str
                if debug:
                    log.warning(f'{key}={value!r},')
                data_dict[key] = value
            return data_dict

        data_dict = _parse(slines, debug=False)
        try:
            qual_info = QualInfo(**data_dict)
        except TypeError as exception:
            qual_info = QualInfo()
            log.debug(f'{db_key: 7d} {qual_str}')
            data_dict = _parse(slines, debug=True)
            log.error(str(exception))
            #raise
        return qual_info
