include Config.mak

# C++ files
CPP_SRC=parmec.cpp input.cpp output.cpp tasksys.cpp mem.cpp map.cpp mesh.cpp timeseries.cpp joints.cpp

# ISPC files
ISPC_SRC=parmec.ispc partition.ispc condet.ispc forces.ispc dynamics.ispc shapes.ispc obstacles.ispc restrain.ispc

# ISPC target
ISPC_TARGET=avx2-i32x8

# Library name
LIB=libparmec

# Program name
EXE=parmec

# Do the rest
include common.mk
