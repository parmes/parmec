# find parmec path
import os, sys
def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None
path = where('parmec4')
if path == None:
  print('ERROR: parmec4 not found in PATH!')
  print('       Download and compile parmec;')
  print('add parmec directory to PATH variable;')
  sys.exit(1)
print('(Found parmec4 at:', path + ')')
sys.path.append(os.path.join (path, 'python'))
from progress_bar import * # and import progress bar
from scipy import spatial  # import scipy
import numpy as np # and numpy

# command line arguments
av = ARGV()
if '-h' in av or '--help' in av:
  print('Beam-like spring-system example:', end=' ')
  print('cantilever beam fixed at x-far-end')
  print('Unit cubes interact via springs')
  print('connected within a radius of influence')
  print('Available arguments:')
  print('  -nx int --> x resolution (or 10)')
  print('  -ny int --> y resolution (or 5)')
  print('  -nz int --> z resolution (or 5)')
  print('  -du float --> duration (or 5.)')
  print('  -st float --> time step (or auto)')
  print('  -ra float --> spring influence radius (or 2.)')
  print('  -h or --help --> print this help')
  sys.exit(0)

# input parameters
nx = int(av[av.index('-nx')+1]) if '-nx' in av else 10
ny = int(av[av.index('-ny')+1]) if '-ny' in av else 5
nz = int(av[av.index('-nz')+1]) if '-nz' in av else 5
du = float(av[av.index('-du')+1]) if '-du' in av else 5.
st = float(av[av.index('-st')+1]) if '-st' in av else -1
ra = float(av[av.index('-ra')+1]) if '-ra' in av else  2.

# materials
matnum = MATERIAL (1E3, 1E9, 0.25)
spring = [-1,-1E6, 1,1E6]
dratio = 10.

# (nx,ny,nz) array of unit cubes
iend = nx*ny*nz-1
progress_bar(0, iend, 'Adding particles:')
x, y, z = np.mgrid[0:nx, 0:ny, 0:nz]
data = zip(x.ravel(), y.ravel(), z.ravel())
datarange = range (0, len(data))
for i in datarange:
  p = data[i]
  nodes = [p[0]-.5, p[1]-.5, p[2]-.5,
	   p[0]+.5, p[1]-.5, p[2]-.5,
	   p[0]+.5, p[1]+.5, p[2]-.5,
	   p[0]-.5, p[1]+.5, p[2]-.5,
	   p[0]-.5, p[1]-.5, p[2]+.5,
	   p[0]+.5, p[1]-.5, p[2]+.5,
	   p[0]+.5, p[1]+.5, p[2]+.5,
	   p[0]-.5, p[1]+.5, p[2]+.5]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  parnum = MESH (nodes, elements, matnum, 0)
  progress_bar(i, iend, 'Adding particles:')

# connecting springs within radius
def add(a,b): return (a[0]+b[0],a[1]+b[1],a[2]+b[2])
def mul(a,b): return (a[0]*b,a[1]*b,a[2]*b)
progress_bar(0, iend, 'Adding springs:')
tree = spatial.KDTree(data)
for i in datarange:
   p = data[i]
   adj = tree.query_ball_point(np.array(p), ra)
   for j in [k for k in adj if k < i]:
     q = data[j]
     x = mul(add(p,q),.5)
     sprnum = SPRING (i, x, j, x, spring, dratio)
   progress_bar(i, iend, 'Adding springs:')

# fixed at x-far-end
for i in datarange[-ny*nz:]:
  RESTRAIN (i, [1,0,0,0,1,0,0,0,1], [1,0,0,0,1,0,0,0,1])

# gravity acceleration
GRAVITY (0., 0., -9.8)

# time step
hc = CRITICAL(perparticle=10)
if st < 0: st = 0.5 * hc[0][0]

# print out statistics
print('%dx%dx%d=%d particles and %d springs' % (nx,ny,nz,parnum,sprnum))
print('10 lowest-step per-particle tuples (critical step, particle index, circular frequency, damping ratio):')
print(hc)
print('Running %d steps of size %g:' % (int(du/st),st))

# run simulation
DEM (du, st, (0.05, 0.01))
