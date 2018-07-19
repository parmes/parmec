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

# Optional paths (comment out {MED/STRUM}INC path if
# a given library is if not available or not needed):
# ===================================================

# MED paths (used as an optional output format)
MEDINC=-I/Users/tomek/Devel/med-3.2.0/build/include
MEDLIB=-L/Users/tomek/Devel/med-3.2.0/build/lib -lmed

# MUMPS path
#MUMPS=/Users/tomek/Devel/MUMPS_5.1.2
MUMPSLIB= -L$(MUMPS)/libseq -lmpiseq -L$(MUMPS)/lib -ldmumps -lmumps_common -lpord -lgfortran
BLASLIB=-L/Users/tomek/Devel/OpenBLAS -lopenblas

# SuiteSparse path
#SUITESPARSE=/Users/tomek/Devel/SuiteSparse

# STRUMPACK paths (used by joints to efficiently solve linear systems)
#STRUMINC=-I/Users/tomek/Devel/STRUMPACK/local/include
STRUMLIB=-L/Users/tomek/Devel/STRUMPACK/local/lib -lstrumpack 
# METIS includes used to compile STRUMPACK
METISINC=-I/Users/tomek/Devel/metis-5.1.0/local/include
# METIS libraries linked with STRUMPACK
METISLIB=-L/Users/tomek/Devel/metis-5.1.0/local/lib -lmetis
# MPI includes used to compile STRUMPACK
MPIINC=$(shell mpicxx --showme:compile)
# MPI libraries linked with STRUMPACK
MPILIB=$(shell mpicxx --showme:link) $(shell mpif90 --showme:link)
# LAPACK libraries linked with STRUCMPACK
LAPACK=-L/opt/local/lib/lapack/ -llapack -lblas -lcblas
# SCALAPACK libraries linked with STRUMPACK
SCALAPACK=-L/opt/local/lib -lscalapack
