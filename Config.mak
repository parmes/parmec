# C++ compiler (ISPC is assumed to be in the PATH; http://ispc.github.io)
CXX=g++

# Debug version
DEBUG=yes

# Non-optional paths:
# ==================

# Python paths (used to parse input files)
PYTHONINC=$(shell python3-config --includes)
PYTHONLIB=$(shell python3-config --ldflags --embed)

# HDF5 paths (used to write output files)
HDF5INC=$(shell pkg-config --cflags hdf5)
HDF5LIB=$(shell pkg-config --libs hdf5) -lhdf5_hl

# Optional paths:
# (comment out if not needed)
# ===========================

# MED paths (used as an optional output format; see also
# http://www.salome-platform.org/user-section/about/med)
#MEDINC=-I/Users/tomek/Devel/med-3.2.0/build/include
#MEDLIB=-L/Users/tomek/Devel/med-3.2.0/build/lib -lmed

# SuiteSparse (when enabled the SuiteSparseQR solver replaces the default, OpenMP-parallel
# skyline LU linear solver, used to calculate joints forces; when your model includes many
# independent systems of particles connected via non-redundant joints, the default solver
# is typically more effective; when your model includes overdetermined (redundant) particle-and-joint
# systems and the default solver fails, then the QR factorisation may be able to provide a solution)
#SUITESPARSE=/Users/tomek/Devel/SuiteSparse
