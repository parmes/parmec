# PARMEC test --> SPRING general nonlinear capability command test

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

spring = [-1, -2E7, -0.05, -1E7, 0, 0, 0.05, 2E7, 1, 4E7]
unload = [-0.05, -1E7, 0, 0, 0.05, 2E7]

# general nonlinear spring with kinematic hardening
sprkin = SPRING (parnum, (0.5, 0.5, 0.5), -1, (0.5, 0.5, 0.5), spring, direction = (-1, 0, 0), unload=unload, ylim=(-1E7, 2E7))

try:
  from scipy.interpolate import interp1d
except:
  print 'ERROR: SciPy interp1d failed to load -->'
  print '       perhaps SciPy needs to be installed'

#vt = [0, 0.0999, 0.1, 0.2999, 0.3, 0.4999]
#vx = [1,      1,  -1,     -1,   1,      1]
#vt = [0, 0.0999, 0.1, 0.1499, 0.15, 0.4999]
#vx = [-1,    -1,   1,      1,   -1,     -1]
# ------------------------------------------------------------------
# loading, unloading, re-loading, loading, unloading, crossing zero, -->
# --> loading, unloading, re-loading, loading, unloading, crossing zero, loading
vt = [0, 0.0999, 0.1, 0.1299, 0.13, 0.2999, 0.3, 0.4999, 0.5, 0.5299, 0.53, 0.7999, 0.8, 1.2]
vx = [1,      1,  -1,     -1,    1,      1,  -1,     -1,   1,      1,   -1,     -1,   1,   1]
fv = interp1d (vt, vx)

def linvel(t): return (float(fv(t)), 0.0, 0.0)

PRESCRIBE (parnum, linvel);

t = HISTORY('TIME') 
dx = HISTORY ('DX', sprkin)
stroke = HISTORY ('STROKE', sprkin)
force = HISTORY ('SF', sprkin)

h = 0.001 #0.9 * CRITICAL()

print 'Time step:', h

DEM (vt[-1], h, (0.05, h))

print 'Generating stroke-force plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, dx, label='dx')
  plt.plot (t, stroke, label='stroke')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('distance $(m)$')
  plt.legend (loc = 'upper right')
  plt.savefig ('tests/spring_general_dxs.png')

  plt.clf ()
  plt.plot (stroke, force)
  plt.xlim ((min(stroke), max(stroke)))
  plt.xlabel ('stroke $(m)$')
  plt.ylabel ('force $(N)$')
  plt.title ('kinematic hardening')
  plt.savefig ('tests/spring_general_sf.png')

except:
  print 't = ', t
  print 'dx = ', dx
  print 'stroke = ', stroke
  print 'force = ', force
