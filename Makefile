# C/C++ compiler
CXX=g++

# C++ files
CPP_SRC=parmec.cpp input.cpp output.cpp tasksys.cpp mem.cpp map.cpp mesh.cpp

# ISPC files
ISPC_SRC=parmec.ispc partition.ispc condet.ispc forces.ispc dynamics.ispc shapes.ispc obstacles.ispc constrain.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# Python paths
PYTHONINC=-I/usr/include/python2.7
PYTHONLIB=-L/usr/lib -lpython2.7

# Program name
EXE=parmec

# Floating point type
REAL=float

# Debug version
DEBUG=yes

# Do the rest
include common.mk
