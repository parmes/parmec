# C++ compiler (ISPC is assumed to be in the PATH; http://ispc.github.io)
CXX=g++

# Python paths
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHONLIB=-L/opt/local/lib -lpython2.7

# HDF5 paths
HDF5INC=-I/opt/local/lib/hdf5-18/include
HDF5LIB=-L/opt/local/lib/hdf5-18/lib -lhdf5 -lhdf5_hl

# MED paths (comment out if not available or needed)
MEDINC=-I/Users/tomek/Devel/med-3.2.0/build/include
MEDLIB=-L/Users/tomek/Devel/med-3.2.0/build/lib -lmed

# Debug version
DEBUG=yes
