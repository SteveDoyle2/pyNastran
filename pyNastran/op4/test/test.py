import pyNastran.op4.op4 as OP4
#print "f = ",op4.__file__


def pass_test1():
    fh = OP4.File(r'C:\Users\steve\Desktop\pyNastran\pyNastran\op4\test\mat_b_dn.op4','r')
    fh.print_header()
    #print fh.nmat = 9
    
    # crash
    a,b,c = fh.Load(nmat=3,skip=0)
    print a
    print b
    print c

def failed_test1():
    fh = OP4.File('mat_b_dn.op4','r')
    fh.print_header()
    #print fh.nmat = 9
    
    a,b,c = fh.Load(nmat=3,skip=0)
    print a
    print b
    print c

def pass_test2():
    fh = OP4.File(r'C:\Users\steve\Desktop\pyNastran\pyNastran\op4\test\mat_b_dn.op4','r')
    #print fh.nmat = 9

    # crash with "unnamed is sparse, skipping for now"
    (a,b,c,d, f,g,h,i) = fh.Load(nmat=9,skip=0) 
    print a
    print b
    print c

def failed_test2():
    fh = OP4.File('mat_b_dn.op4','r')
    #print fh.nmat = 9

    # ValueError:  need more than 8 values to unpack
    (a,b,c,d, f,g,h,i) = fh.Load(nmat=9,skip=0) 
    print a
    print b
    print c

print "-------------"
failed_test1()
print "*********"
failed_test2()

pass_test1()
pass_test2()
