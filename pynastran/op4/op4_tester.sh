#!/bin/bash

export PYTHONPATH="${PYTHONPATH}:/my/other/path"

export PYTHONPATH="${PYTHONPATH}:/home/git/pysa/pynastran/pynastran/op4:/home/git/pysa/binaryornot/binaryornot"
export PYTHON="python3"
#set PYFILEPATH=""
#%PYTHON% "op4.py"
#%PYTHON% "op4_bis.py"
#%PYTHON% -i "C:\Users\ecannizzaro1\Documents\git\pysa\pynastran\pynastran\op4\op4_temp.py"
#%PYTHON% -i "C:\Users\ecannizzaro1\Documents\git\pysa\pynastran\pynastran\op4\op4_utils.py"
${PYTHON} op4_tester.py
#%PYTHON% -i "C:\Users\ecannizzaro1\Documents\git\pysa\pynastran\pynastran\op4\test_op4.py"
#%PYTHON% -i "C:\Users\ecannizzaro1\Documents\git\pysa\pynastran\pynastran\op4\demo_op4.py"
#%PYTHON% -i "C:\Users\ecannizzaro1\Documents\git\pysa\pynastran\pynastran\op4\demo_op4_all_files.py"
