# PARMEC test --> ANALYTICAL particle

matnum = MATERIAL (1E3, 1E9, 0.25)

parnum = ANALYTICAL ([1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 1.0,
  [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], (0.0, 0.0, 0.0), matnum)

VELOCITY (parnum, angular = (0, 0, 1))

DEM (10, 0.01, 0.1)
