"""
Note: In some OFP table descriptions (OEE, OEF, OES for example), you
will see formats such as ACODE,4=05, or TCODE,1=02
(versus ACODE=05 or TCODE=02). The integer values 4 and 1 in these
examples are function codes. Function codes specify operations to
perform on the value in the data block. The operation result will
then be used to determine the data format.

The following lists the available function codes and their operation:

Function Codes Operation
============== =========
  1            if (item_name/1000 = 2,3,6) then return 2, else return 1
  2            mod(item_name,100)
  3            mod(item_name,1000)
  4            item_name/10
  5            mod(item_name,10)
  6            if iand(item_name,8) != then set to 0, else set to 1
  7            if item_name/1000
               = 0 or 2, then set to 0
               = 1 or 3, then set to 1
               > 3, then set to 2.
  >65535       iand(item_name,iand(func_code,65535))

For example, if a value of 100 is found in an ACODE,4 field, the
function_code of 4 results in the operation:
    item_name/10 = 100/10 = 10.

Thus the data format under the ACODE,4=10 row would be used.
"""
import warnings

def func1(item_code):
    if item_code // 1000 in [2, 3, 6]:
        return 2
    return 1

def func2(item_code):
    return item_code % 100

def func3(item_code):
    return item_code % 1000

def func4(item_code):
    return item_code // 10

def func5(item_code):
    return item_code % 10

def func6(item_code):
    warnings.warn('Function code 6 method not verified',
                  RuntimeWarning)
    if item_code & 8:
        return 0
    return 1

def func7(item_code):
    v = item_code // 1000
    if v in [0, 2]:
        return 0
    if v in [1, 3]:
        return 1
    return 2

def funcbig(func_code, item_code):
    return item_code & (func_code & 65535)

#self._code_funcs = {
    #1: func1, 2: func2, 3: func3, 4: func4,
    #5: func5, 6: func6, 7: func7,
    #'big': funcbig,
#}

def func7(item_code: int) -> int:
    v = item_code // 1000
    if v in [0, 2]:
        return 0
    if v in [1, 3]:
        return 1
    return 2
