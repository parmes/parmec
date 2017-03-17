print '---------------------------------'
print 'Parmec acceleration sweep example'
print '---------------------------------'

# input parameters
N = 10 # number of bodies along one direction
gap = 0.001 # gap between cubes
lofq = 1 # low frequency
hifq = 5 # high frequency
amag = 1 # acceleration magnitude
stop = 10 # simulation duration
step = 1e-04 # time step estimation
ratio = 0.2 # critical time step ratio

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

# prepare acceleration sweep based input motion
from acc_sweep import *
(vt, vd, vv, va) = acc_sweep (step, stop, lofq, hifq, amag)
tsv = [None]*(len(vt)+len(vd))
tsv[::2] = vt
tsv[1::2] = vv
tsv = TSERIES (tsv)
ts0 = TSERIES (0.0)
linvel = (tsv, tsv, tsv)
angvel = (ts0, ts0, ts0)

# cube material
matnum = MATERIAL (100, 1E6, 0.25)

# cube creation function
def cube (x, y, z):
  nodes = [x+0.0, y+0.0, z+0.0,
	   x+0.1, y+0.0, z+0.0,
	   x+0.1, y+0.1, z+0.0,
	   x+0.0, y+0.1, z+0.0,
	   x+0.0, y+0.0, z+0.1,
	   x+0.1, y+0.0, z+0.1,
	   x+0.1, y+0.1, z+0.1,
	   x+0.0, y+0.1, z+0.1]
  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
  parnum = MESH (nodes, elements, matnum, colors)
  CONSTRAIN (parnum, angular=[1, 0, 0, 0, 1, 0, 0, 0, 1])
  ANALYTICAL (particle=parnum)
  return parnum

# generate NxNxN array of cubes
parnum = 0
ijkmap = {}
for i in range (0,N):
  for j in range (0,N):
    for k in range (0,N):
      parnum = cube (i*(0.1+gap), j*(0.1+gap), k*(0.1+gap))
      ijkmap[(i,j,k)] = parnum

# prescribe outer shell motion
for (i,j,k) in ijkmap:
  outer = [0, N-1]
  if i in outer or j in outer or k in outer:
    num = ijkmap[(i,j,k)]
    PRESCRIBE (num, linvel, angvel) # outer most shell of bodies

# define contact spring curves
spring_curve = [-1-gap, -1E3, -gap, 0, 1, 0]
damper_curve = [-1, -7, 1, 7]

# insert contact springs
sprnum = 0
ijkmax = N-1
for (i,j,k) in ijkmap:
  if i < ijkmax:
    p1 = (i*(0.1+gap)+0.1, j*(0.1+gap)+0.05, k*(0.1+gap)+0.05)
    p2 = (i*(0.1+gap)+0.1+gap, j*(0.1+gap)+0.05, k*(0.1+gap)+0.05)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i+1,j,k)]
    sprnum = SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (1, 0, 0))
  if j < ijkmax:
    p1 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1, k*(0.1+gap)+0.05)
    p2 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.1+gap, k*(0.1+gap)+0.05)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i,j+1,k)]
    sprnum = SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (0, 1, 0))
  if k < ijkmax:
    p1 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.05, k*(0.1+gap)+0.1)
    p2 = (i*(0.1+gap)+0.05, j*(0.1+gap)+0.05, k*(0.1+gap)+0.1+gap)
    n1 = ijkmap[(i,j,k)]
    n2 = ijkmap[(i,j,k+1)]
    sprnum = SPRING (n1, p1, n2, p2, spring_curve, damper_curve, (0, 0, 1))

# set gravity
GRAVITY (0.0, 0.0, -9.81)

# request time histories
ths = HISTORY ('TIME')
vxh = HISTORY ('VX', ijkmap[(N/2,N/2,N/2)])
vyh = HISTORY ('VY', ijkmap[(N/2,N/2,N/2)])
vzh = HISTORY ('VZ', ijkmap[(N/2,N/2,N/2)])
dxh = HISTORY ('DX', ijkmap[(N/2,N/2,N/2)])
dyh = HISTORY ('DY', ijkmap[(N/2,N/2,N/2)])
dzh = HISTORY ('DZ', ijkmap[(N/2,N/2,N/2)])

# print statistics
print parnum, 'bodies'
print sprnum, 'springs'

# estimate critical time step
hcrit = CRITICAL()
step = ratio * hcrit
print 'Critical time step: %.2e' % hcrit

# run simulation
print 'Running', int(stop/step), 'DEM steps of size %.2e' % step, '...'
t = DEM (stop, step, (0.1, 0.01))
print 'Finished after %f seconds' % t

# generate plots
print 'Generating time history plots...'
vhs = []
dhs = []
for i in range(0, len(ths)):
  vhs.append((vxh[i]**2+vyh[i]**2+vzh[i]**2)**0.5)
  dhs.append((dxh[i]**2+dyh[i]**2+dzh[i]**2)**0.5)

try:
  import matplotlib.pyplot as plt
  plt.clf ()
  plt.plot (ths, vhs)
  plt.xlim ((0, ths[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('|velocity| $(m/s)$')
  plt.savefig ('examples/asweep_velo.png')

  plt.clf ()
  plt.plot (ths, dhs)
  plt.xlim ((0, ths[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('|displacement| $(m)$')
  plt.savefig ('examples/asweep_disp.png')

except:
  print 'time = ', ths
  print 'velo = ', vhs
  print 'disp = ', dhs
