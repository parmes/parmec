# PARMEC test --> DEM (...,interval,...)

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

VELOCITY (parnum, angular = (0, 0, 1))

t = HISTORY ('TIME')
o = HISTORY ('OZ', parnum)

try:
  from scipy.interpolate import interp1d
except:
  print 'ERROR: SciPy interp1d failed to load -->'
  print '       perhaps SciPy needs to be installed'

dt_files = interp1d([0,5,10], [0.5,3,0.2])
dt_hist = interp1d([0,5,10], [2,0.1,1])

DEM (10, 0.1, interval=(dt_files, dt_hist))

print 'Generating plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, o, 'rs')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('$\omega_z$ $(m)$')
  plt.title ('angular velocity')
  plt.savefig ('tests/interval.png')

except:
  print 't = ', t
  print 'o = ', o
