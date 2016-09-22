# PARMEC test --> PRESCRIBE command test

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

def linvel(t): return (t, 0, 0)
def linacc(t): return (1, 0, 0)
def angvel(t): return (0, 0, t)
def angacc(t): return (0, 0, 1)

PRESCRIBE (parnum, linvel, angvel)
# PRESCRIBE (parnum, linacc, angacc, 'aa')
# PRESCRIBE (parnum, linvel, angacc, 'va')
# PRESCRIBE (parnum, linacc, angvel, 'av')

t = HISTORY ('TIME')
vx = HISTORY ('VX', parnum)
oz = HISTORY ('OZ', parnum)

DEM (2.0, 0.01, (0.1, 0.01))

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, vx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('vx $(m/s)$')
  plt.savefig ('tests/prescribe_vx.png')

  plt.clf ()
  plt.plot (t, oz)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('oz $(rad/s)$')
  plt.savefig ('tests/prescribe_oz.png')

except:
  print 'time = ', t
  print 'vx = ', vx
  print 'oz = ', oz
