from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
from numpy import matrix


def printBadPath(path):
    msg = ''
    fullpath = os.path.abspath(path)
    while fullpath:
        if os.path.exists(fullpath):
            msg += 'passed:  %s\n' % (fullpath)
        else:
            msg += 'failed:  %s\n' % (fullpath)
        fullpath = os.path.dirname(fullpath)
    return msg


def get_files_of_type(dirname, extension='.txt', maxSize=100.):
    """
    gets all the files in the specified directory with a given extension
    @param dirname the directory name
    @param extension list of filetypes to get (default='.txt')
    @param maxSize size in MB for max file size
    """
    files = os.listdir(dirname)
    files2 = []
    for fname in files:
        (f, ext) = os.path.splitext(fname)
        if extension in ext:
            f = os.path.join(dirname, fname)
            if os.path.getsize(f) / (1024. * 1024.) <= maxSize:  # convert to MB
                files2.append(f)
    return files2


def deObscure(num):
    """
    unpacks an "obscured" number...similar to binary, but base 52
    A base 52 value takes up a fewer characters than a base 10 number
    which helps to do Mat12345678 when there's an 8 character limit on variable names.
    """
    print("***")
    print("type(num) = ", type(num))
    num.reverse()
    vals = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ]
    #vals = ['0','1']
    dictA = {}
    n = len(vals)
    for i in xrange(n):
        dictA[vals[i]] = i

    print("n = ", n)
    val = 0
    for i, letter in enumerate(reversed(num)):
        print("letter = ", letter)
        val += dictA[letter] * n ** i
        print("factor = ", dictA[letter] * n ** i)
    print("val = ", val)
    return val


def obscure(num):
    """
    takes a large number and shrinks it down...similar to binary, but base 52
    A base 52 value takes up a fewer characters than a base 10 number
    which helps to do Mat12345678 when there's an 8 character limit on variable names.
    """
    lenNum = len(str(num))
    vals = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', ]
    #vals = ['0','1']
            #'0','1','2','3','4','5','6','7','8','9']  # 0,1,2,...
    dictA = {}
    n = len(vals)
    for i in xrange(n):
        dictA[i] = vals[i]

    pack = []
    i = 0
    while num > 0:
        print("num = %s" % (num))
        print("factor = %s" % (num % n))
        var = dictA[num % n]
        num = num / n
        pack.append(var)
        i += 1
        if i == 100:
            break
    print(pack)
    print('')
    print("%s > %s" % (lenNum, len(pack)))
    return pack


def is_binary(filename):
    """
    Return true if the given filename is binary.
    @raise EnvironmentError: if the file does not exist or cannot be accessed.
    @attention: found @ http://bytes.com/topic/python/answers/21222-determine-file-type-binary-text on 6/08/2010
    @author: Trent Mick <TrentM@ActiveState.com>
    @author: Jorge Orpinel <jorge@orpinel.com>
    @warning this may not work for unicode...
    """
    fin = open(filename, 'rb')
    try:
        CHUNKSIZE = 1024
        while 1:
            chunk = fin.read(CHUNKSIZE)
            if b'\0' in chunk:  # found null byte
                return True
            if len(chunk) < CHUNKSIZE:
                break  # done
    finally:  # no need for ???
        fin.close()

    return False


def ListPrint(listA):
    """
    Prints a list, numpy array, or numpy matrix in an abbreviated format.
    Useful for debugging.
    @param listA list, numpy array, or numpy matrix
    @retval msg the clean string representation of the object
    """
    if len(listA) == 0:
        return '[]'

    msg = '['
    if isinstance(listA, matrix):
        (nrows, ncols) = listA.shape
        for irow in xrange(nrows):
            msg += '['
            for icol in xrange(ncols):
                msg += '%-10g,' % (listA[irow, icol])
            msg = msg[:-1]
            msg += '],\n '
        msg = msg[:-1]
        msg += ']'
    else:
        for a in listA:
            #print "a = ", a, type(a)
            if isinstance(a, unicode):
                msg += ' %s,' % (str(a))
            elif a is None:
                msg += ' None,'
            elif isinstance(a, float):
                msg += ' %-4.2f,' % (a)
            elif isinstance(a, int):
                msg += ' %g,' % (a)
            else:
                try:
                    msg += ' %g,' % (a)
                except TypeError:
                    print("a = |%s|" % (a))
                    raise
                ###
            ###
        ###
        msg = msg[:-1]
        msg += ' ]'
    ###
    return msg

if __name__ == '__main__':
    n = 99999999
    o = obscure(n)
    print(''.join(o))
    deObscure(o)
