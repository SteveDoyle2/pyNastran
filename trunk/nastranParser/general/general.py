
def ListPrint(listA):
    if len(listA)==0:
        return '[]'
    ###

    msg = '['
    for a in listA:
        if isinstance(a,str):
            msg += ' %s,' %(a)
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
    return msg

