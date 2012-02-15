from numpy import matrix
from types import NoneType

def deObscure(num):
    """
    unpacks an "obscured" number...similar to binary, but base 52
    """
    print "***"
    print "type(num) = ",type(num)
    num.reverse()
    vals = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
            'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',]
    #vals = ['0','1']
    dictA = {}
    n = len(vals)
    for i in range(n):
        dictA[vals[i]] = i

    print "n = ",n
    val = 0
    for i,letter in enumerate(reversed(num)):
        print "letter = ",letter
        val += dictA[letter]*n**i
        print "factor = ",dictA[letter]*n**i
    print "val = ",val
    return val

def obscure(num):
    """
    takes a large number and shrinks it down...similar to binary, but base 52
    """
    lenNum = len(str(num))
    vals = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',
            'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',]
    #vals = ['0','1']
            #'0','1','2','3','4','5','6','7','8','9']  # 0,1,2,...
    dictA = {}
    n = len(vals)
    for i in range(n):
        dictA[i] = vals[i]

    pack = []
    i=0
    while num>0:
        print "num = ",num
        print "factor = ",num%n
        var = dictA[num%n]
        num = num/n
        pack.append(var)
        i+=1
        if i==100:
            break
    print pack
    print 
    print "%s > %s" %(lenNum,len(pack))
    return pack
    

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

if __name__=='__main__':
    n = 99999999
    o = obscure(n)
    print ''.join(o)
    deObscure(o)