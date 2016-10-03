# PARMEC test --> PRESCRIBE command test (velocity sweep)
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

stop = 5.0

def linvel(t):
  amag = 1.0
  lofq = 0.0
  hifq = 5.0
  # derivative of a = amag * sin (2.0*pi*(lofq+(hifq-lofq)*t/stop)*t) -->
  v = amag * cos (2.0*pi*(lofq+(hifq-lofq)*t/stop)*t) * (2.0*pi*lofq + 4.0*pi*(hifq-lofq)*t/stop)
  return (v, 0, 0)

PRESCRIBE (parnum, linear = linvel)

t = HISTORY ('TIME')
vx1 = HISTORY ('VX', parnum)
dx = HISTORY ('DX', parnum)

DEM (stop, 0.001, (0.05, 0.01))

vx0 = []
for s in t: vx0.append(linvel(s)[0])

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, vx0, label = 'input', linestyle = '--', marker = '.')
  plt.plot (t, vx1, label = 'output')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('vx $(m/s)$')
  plt.legend()
  plt.savefig ('tests/prescribe_vsweep_vx.png')

  plt.clf ()
  plt.plot (t, dx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dx $(m)$')
  plt.savefig ('tests/prescribe_vsweep_dx.png')

except:
  print 'time = ', t
  print 'vx0 = ', vx0
  print 'vx1 = ', vx1
  print 'dx = ', dx
