import os, sys

def where(program):
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, program)):
      return path
  return None

path = where('parmec4')

if path == None:
  print 'ERROR: parmec4 not found in PATH!'
  sys.exit(1)

print '(Found parmec4 at:', path + ')'

print '=================================='
print 'PARMEC: acceleration sweep example'
print '=================================='

sys.path.append(os.path.join (path, 'python'))

from acc_sweep import *

stop = 5.0   # duration of the simulation
ostp = 0.05  # output step
lofq = 5     # low frequency for the sweep
hifq = 15    # high frequency for the sweep
amag = 10.0  # acceleration magnitude
nbod = 2     # number of bodies
nele = 2     # number of elements per body (along x, y, z)
l = 0.1      # length of one body
w = 0.1      # widhth of one body
h = 0.1      # height of one body
gap = 0.002  # gap

def nodes(x,y,z): # cube corner nodes
  return [x, y, z,
	  x+w, y, z,
	  x+w, y+l, z,
	  x, y+l, z,
	  x, y, z+h,
	  x+w, y, z+h,
	  x+w, y+l, z+h,
	  x, y+l, z+h]

matnum = MATERIAL (1E3, 1E9, 0.25)

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

parnum = []
for i in range (0, nbod+2):
  j = MESH (nodes(0, i*(l+gap), 0), elements, matnum, i)
  parnum.append (j)

spring = [-1.0-gap, -1E5, -gap, 0.0, 1.0, 0.0]
damper = [-1.0, -0.9E3, 1.0, 0.9E3]

for i in range (0, nbod+1):
  p0 =  (0, i*(l+gap)+l, 0)
  p1 =  (0, (i+1)*(l+gap), 0)
  SPRING (parnum[i], p0, parnum[i+1], p1, spring, damper, (0, 1, 0))
  p0 =  (w, i*(l+gap)+l, 0)
  p1 =  (w, (i+1)*(l+gap), 0)
  SPRING (parnum[i], p0, parnum[i+1], p1, spring, damper, (0, 1, 0))
  p0 =  (w, i*(l+gap)+l, h)
  p1 =  (w, (i+1)*(l+gap), h)
  SPRING (parnum[i], p0, parnum[i+1], p1, spring, damper, (0, 1, 0))
  p0 =  (0, i*(l+gap)+l, h)
  p1 =  (0, (i+1)*(l+gap), h)
  SPRING (parnum[i], p0, parnum[i+1], p1, spring, damper, (0, 1, 0))

step = 0.2 * CRITICAL()

print 'Time step:', step

(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)

try:
  from scipy.interpolate import interp1d
except:
  print 'ERROR: SciPy interp1d failed to load -->'
  print '       perhaps SciPy needs to be installed'
  sys.exit(1)

vel = interp1d(vt, vv) # acceleration sweep linear spline

def linvel(t): # acceleration sweep function
  return (0, vel(t), 0)

def angvel(t): return (0, 0, 0) # prescribed angular velocity

PRESCRIBE (parnum[0], linvel, angvel, kind = 'vv')
PRESCRIBE (parnum[nbod+1], linvel, angvel, kind = 'vv')

t = HISTORY ('TIME')
vx = HISTORY ('VX', parnum[1]) # FIXME: these time histories are zero upon termination
dx = HISTORY ('DX', parnum[1])

DEM (stop, step, ostp)

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, vx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('vx $(m/s)$')
  plt.savefig ('examples/asweep_vx.png')

  plt.clf ()
  plt.plot (t, dx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dx $(m)$')
  plt.savefig ('examples/asweep_dx.png')

except:
  print 'time = ', t
  print 'vx = ', vx
  print 'dx = ', dx
