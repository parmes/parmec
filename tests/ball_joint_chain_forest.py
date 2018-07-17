# PARMEC test --> BALL_JOINT command test with many indepedent chains
#                 testing parallel performance of the linear system solve

nedge  = 10
nchain = 10

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

for x in range (0, nchain*nedge, nchain):
  for y in range (0, nchain*nedge, nchain):
    for i in range (0, nchain):
      nodes = [x+i, y+i, i,
	       x+i+1, y+i, i,
	       x+i+1, y+i+1, i,
	       x+i, y+i+1, i,
	       x+i, y+i, i+1,
	       x+i+1, y+i, i+1,
	       x+i+1, y+i+1, i+1,
	       x+i, y+i+1, i+1]

      elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

      colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

      parnum = MESH (nodes, elements, matnum, colors)

      if i: BALL_JOINT (parnum, (x+i, y+i, i), prevpar)

      prevpar = parnum

    sprnum = BALL_JOINT (parnum, (x+i+1, y+i+1, i+1))

GRAVITY (0., 0., -10.)

h = 0.001
DEM (5.0, h, (0.05, h))
