# PARMEC test --> CONSTRAIN command test

matnum = MATERIAL (1E3, 1E9, 0.25)

for i in range (0,7):
  nodes = [2*i, 2*i, 0,
	   2*i+1, 2*i, 0,
	   2*i+1, 2*i+1, 0,
	   2*i, 2*i+1, 0,
	   2*i, 2*i, 1,
	   2*i+1, 2*i, 1,
	   2*i+1, 2*i+1, 1,
	   2*i, 2*i+1, 1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

  parnum = MESH (nodes, elements, matnum, i)

  VELOCITY (parnum, (1, 1, 1), (1, 1, 1))

  if i == 1:
    CONSTRAIN (parnum, [1, 0, 0])
  if i == 2:
    CONSTRAIN (parnum, [1, 0, 0, 0, 1, 0])
  if i == 3:
    CONSTRAIN (parnum, [1, 0, 0, 0, 1, 0, 0, 0, 1])
  if i == 4:
    CONSTRAIN (parnum, angular = [1, 0, 0])
  if i == 5:
    CONSTRAIN (parnum, angular = [1, 0, 0, 0, 1, 0])
  if i == 6:
    CONSTRAIN (parnum, angular = [1, 0, 0, 0, 1, 0, 0 ,0, 1])

DEM (5, 0.1)
