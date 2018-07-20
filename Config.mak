# C++ compiler (ISPC is assumed to be in the PATH; http://ispc.github.io)
CXX=g++

# Debug version
DEBUG=yes

# Non-optional paths:
# ==================

# Python paths (used to parse input files)
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHONLIB=-L/opt/local/lib -lpython2.7

# HDF5 paths (used to write output files)
HDF5INC=-I/opt/local/lib/hdf5-18/include
HDF5LIB=-L/opt/local/lib/hdf5-18/lib -lhdf5 -lhdf5_hl

# Optional paths:
# (comment out if not needed)
# ===========================

# MED paths (used as an optional output format; see also
# http://www.salome-platform.org/user-section/about/med)
MEDINC=-I/Users/tomek/Devel/med-3.2.0/build/include
MEDLIB=-L/Users/tomek/Devel/med-3.2.0/build/lib -lmed

# SuiteSparse (when enabled the SuiteSparseQR solver replaces the default, OpenMP-parallel
# skyline LU linear solver, used to calculate joints forces; when your model includes many
# independent systems of particles connected via non-redundant joints, the default solver
# is typically more effective; when your model includes overdetermined (redundant) particle-and-joint
# systems and the default solver fails, then the QR factorisation may be able to provide a solution)
#SUITESPARSE=/Users/tomek/Devel/SuiteSparse
