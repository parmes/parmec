# PARMEC test --> global DAMPING using a projectile example
print 'Projectile test...'

v0 = (2.54, 0, 12.7)
w0 = (0, 0, 0)
gravity = 9.81456;
stop = 1.976
step = stop / 1024.0

num = ANALYTICAL ([1, 1, 1, 0, 0, 0], 0.45359237, position=(0, 0, 0))
VELOCITY (num, v0, w0)

def lindamp(t): return (1, 1, 1)
def angdamp(t): return (0, 0, 0)
DAMPING (lindamp, angdamp)

GRAVITY (0., 0., -gravity)

T = HISTORY('TIME')
X = HISTORY('PX', num)
Y = HISTORY('PY', num)
Z = HISTORY('PZ', num)

print 'Calculating...'
DEM (stop, step, (64*step, step))

for i in range(0, len(X)): X[i] = X[i]/0.0254 # conver to inches
for i in range(0, len(Y)): Y[i] = Y[i]/0.0254
for i in range(0, len(Z)): Z[i] = Z[i]/0.0254

print 'Correctness test...',
value = X[-1]
exact = 86.138
error = abs (value - exact) / exact
if (error < 0.005): print 'PASSED'
else:
  print 'FAILED'
  print '(', 'Computed x distance was %.3f' % value, 'while the reference value is %.3f' % exact, ')'

try:
  import matplotlib.pyplot as plt
  print 'Generating trajectory plot...'
  plt.clf ()
  plt.plot (T, X, label='x')
  plt.plot (T, Y, label='y')
  plt.plot (T, Z, label='z')
  plt.axis (xmin = 0, xmax = stop, ymin = 0, ymax = 180)
  plt.xlabel ('Time [s]')
  plt.ylabel ('Position [in]')
  plt.legend (loc = 'upper right')
  plt.savefig ('tests/damping.png')
except ImportError:
  print 'Matplot lib was not installed --> skipping plots'
