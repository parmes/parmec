# looking for the parmec path
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  print '       Download and compile parmec; add parmec directory to PATH variable;'
  sys.exit(1)
print '(Found parmec4 at:', path + ')'
sys.path.append(os.path.join (path, 'python'))
from progress_bar import *
from scipy import spatial
import numpy as np

nx = 10
ny = 10
nz = 10
radius = 1.

# materials
matnum = MATERIAL (1E3, 1E9, 0.25)
spring = [-1,-1E7, 1,1E7]
dashpot = [-1, -8E5, 1, 8E5]

# (nx,ny,nz) array of unit cubes
iend = nx*ny*nz-1
printProgressBar(0, iend, 'Adding particles:')
x, y, z = np.mgrid[0:nx, 0:ny, 0:nz]
data = zip(x.ravel(), y.ravel(), z.ravel())
datarange = range (0, len(data))
for i in datarange:
  p = data[i]
  nodes = [p[0], p[1], p[2],
	   p[0]+1, p[1], p[2],
	   p[0]+1, p[1]+1, p[2],
	   p[0], p[1]+1, p[2],
	   p[0], p[1], p[2]+1,
	   p[0]+1, p[1], p[2]+1,
	   p[0]+1, p[1]+1, p[2]+1,
	   p[0], p[1]+1, p[2]+1]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  parnum = MESH (nodes, elements, matnum, 0)
  printProgressBar(i, iend, 'Adding particles:')

# randomized springs within radius
printProgressBar(0, iend, 'Adding springs:')
tree = spatial.KDTree(data)
for i in datarange:
   adj = tree.query_ball_point(np.array(data[i]), radius)
   for j in adj:
     p = tuple((np.array(data[i])+ np.random.rand(1,3)).ravel())
     q = tuple((np.array(data[j])+ np.random.rand(1,3)).ravel())
     SPRING (i, p, j, q, spring, dashpot)
   printProgressBar(i, iend, 'Adding springs:')

# fix at far end
for i in datarange[-ny*nz:]:
  SPRING (i, data[i], -1, data[i], spring, dashpot)

GRAVITY (0., 0., -10.)

step = 0.05 * CRITICAL()
print 'Time step:', step
DEM (5.0, step, (0.05, step))
