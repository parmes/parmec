# C/C++ compiler
CXX=g++

# C++ files
CPP_SRC=parmec.cpp input.cpp output.cpp tasksys.cpp mem.cpp map.cpp mesh.cpp

# ISPC files
ISPC_SRC=parmec.ispc partition.ispc condet.ispc forces.ispc dynamics.ispc shapes.ispc obstacles.ispc constrain.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# Python paths
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHONLIB=-L/opt/local/lib -lpython2.7

# Debug version
DEBUG=yes

# Library name
LIB=libparmec

# Program name
EXE=parmec

# Do the rest
include common.mk
