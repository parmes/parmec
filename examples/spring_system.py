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

nx = 20
ny = 5
nz = 5
radius = 1.
dura = 0.001

# materials
matnum = MATERIAL (1E3, 1E9, 0.25)
spring = [-1,-1E7, 1,1E7]
dashpot = [-1, -8E5, 1, 8E5]

# (nx,ny,nz) array of unit cubes
iend = nx*ny*nz-1
progress_bar(0, iend, 'Adding particles:')
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
  progress_bar(i, iend, 'Adding particles:')

# randomized springs within radius
def add(a,b): return (a[0]+b[0],a[1]+b[2],a[2]+b[2])
def sub(a,b): return (a[0]-b[0],a[1]-b[2],a[2]-b[2])
def mul(a,b): return (a[0]*b,a[1]*b,a[2]*b)
progress_bar(0, iend, 'Adding springs:')
tree = spatial.KDTree(data)
o5 = (.5,.5,.5)
for i in datarange:
   adj = tree.query_ball_point(np.array(data[i]), radius)
   for j in adj:
     p = add(data[i],o5)
     q = add(data[j],o5)
     x = mul(add(p,q),.5)
     sprnum = SPRING (i, x, j, x, spring, dashpot)
   progress_bar(i, iend, 'Adding springs:')

# fix at far end
for i in datarange[-ny*nz:]:
  RESTRAIN (i, [1,0,0,0,1,0,0,0,1], [1,0,0,0,1,0,0,0,1])

GRAVITY (0., 0., -10.)

step = 0.001 * CRITICAL()
print '%d particles and %d springs' % (parnum, sprnum)
print 'Running %d steps of size %g:' % (int(dura/step), step)
DEM (dura, step, (0.05, 0.01))
