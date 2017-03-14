# C++ compiler (ISPC is assumed to be in the PATH; http://ispc.github.io)
CXX=g++

# C++ files
CPP_SRC=parmec.cpp input.cpp output.cpp tasksys.cpp mem.cpp map.cpp mesh.cpp timeseries.cpp

# ISPC files
ISPC_SRC=parmec.ispc partition.ispc condet.ispc forces.ispc dynamics.ispc shapes.ispc obstacles.ispc constrain.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# Python paths
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHONLIB=-L/opt/local/lib -lpython2.7

# HDF5 paths
HDF5INC=-I/opt/local/include
HDF5LIB=-L/opt/local/lib -lhdf5 -lhdf5_hl

# Debug version
DEBUG=yes

# Library name
LIB=libparmec

# Program name
EXE=parmec

# Do the rest
include common.mk
