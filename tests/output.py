# PARMEC test --> OUTPUT command test

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

for i in range (0, 10):
  nodes = [i, i, i,
	   i+1, i, i,
	   i+1, i+1, i,
	   i, i+1, i,
	   i, i, i+1,
	   i+1, i, i+1,
	   i+1, i+1, i+1,
	   i, i+1, i+1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  parnum = MESH (nodes, elements, matnum, colors)

  if i: SPRING (parnum, (i, i, i), prevpar, (i, i, i), [-1,-1E7, 1,1E7], [-1, -8E5, 1, 8E5])

  prevpar = parnum

SPRING (parnum, (i+1, i+1, i+1), -1, (i+1, i+1, i+1), [-1,-1E7, 1,1E7], [-1, -8E5, 1, 8E5])

GRAVITY (0, 0, -10)
                                      # default output files --> ./tests/output0.vtk.*
OUTPUT (['NUMBER', 'LINVEL'], [0, 1]) # output files --> ./tests/output1.vtk.*
OUTPUT (['NUMBER', 'FORCE', 'TORQUE'], [5, 6, 7]) # output files --> ./tests/output2.vtk.*

h = 0.2 * CRITICAL()

DEM (2.5, h, 0.05)
