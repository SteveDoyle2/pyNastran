g++ grid.cpp coord.cpp rods.cpp shell.cpp main.cpp -o run.exe
swig -python -c++ bdf.i
g++ bdf_wrap.cxx -I F:\Anaconda\include
