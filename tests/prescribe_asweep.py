# PARMEC test --> PRESCRIBE command test (acceleration sweep)
from math import sin, cos, pi

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes = [0, 0, 0,
         1, 0, 0,
	 1, 1, 0,
	 0, 1, 0,
	 0, 0, 1,
	 1, 0, 1,
	 1, 1, 1,
	 0, 1, 1]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

parnum = MESH (nodes, elements, matnum, colors)

def linacc(t):
  amag = 1.0
  lofq = 1.0
  hifq = 10.0
  stop = 10.0
  a = amag * sin (2.0*pi*(lofq+(hifq-lofq)*t/stop)*t)
  return (a, 0, 0)

PRESCRIBE (parnum, linear = linacc, kind = 'av') # FIXME: one directional travel instead of sweep

t = HISTORY ('TIME')
vx = HISTORY ('VX', parnum)
dx = HISTORY ('DX', parnum)

DEM (10.0, 0.01, (0.05, 0.01))

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, vx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('vx $(m/s)$')
  plt.savefig ('tests/prescribe_asweep_vx.png')

  plt.clf ()
  plt.plot (t, dx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dx $(m)$')
  plt.savefig ('tests/prescribe_asweep_dx.png')

except:
  print 'time = ', t
  print 'vx = ', vx
  print 'dx = ', dx
