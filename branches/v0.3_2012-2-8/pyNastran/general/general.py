from numpy import matrix
from types import NoneType

def ListPrint(listA):
    if len(listA)==0:
        return '[]'
    ###

    msg = '['
    if isinstance(listA,matrix):
        (nrows,ncols) = listA.shape
        for irow in range(nrows):
            msg += '['
            for icol in range(ncols):
                msg += '%-10g,' %(listA[irow,icol])
            ###
            msg = msg[:-1]
            msg += '],\n '
        ###
        msg = msg[:-1]
        msg += ']'

    else:
        for a in listA:
            #print "a = ",a,type(a)
            if isinstance(a,str):
                msg += ' %s,' %(a)
            elif isinstance(a,NoneType):
                msg += ' None,'
            elif isinstance(a,float):
                msg += ' %-4.2f,' %(a)
            elif isinstance(a,int):
                msg += ' %g,' %(a)
            else:
                try:
                    msg += ' %g,' %(a)
                except TypeError:
                    print "a = |%s|" %(a)
                    raise
                ###
            ###
        ###
        msg = msg[:-1]
        msg += ' ]'
    ###
    return msg

