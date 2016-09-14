# PARMEC test --> MESH turned into an ANALYTICAL particle
import sys
sys.path.append('python')
from mesh_hex import MESH_HEX

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes = [[0, 0, 0],
         [1, 0, 0],
	 [1, 1, 0],
	 [0, 1, 0],
	 [0, 0, 1],
	 [1, 0, 1],
	 [1, 1, 1],
	 [0, 1, 1]]

colors = [0, 1, 2, 3, 4, 5]

parnum = MESH_HEX (nodes, 5, 5, 5, matnum, colors)

ANALYTICAL (particle = parnum)

VELOCITY (parnum, angular = (0, 0, 1))

DEM (10, 0.1)
