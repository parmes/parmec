# PARMEC test --> SPRING curve correctness tests

MATERIAL (1, 1, 0.25)

parnum = ANALYTICAL ()

# TODO: generate random spring within [-1, -1, 1, 1] box

sprnum = SPRING (parnum, (0, 0, 0), -1, (0, 0, 0), [-1,-1, 1, 1], direction = (-1, 0, 0))

h = 0.02

x = HISTORY ('STROKE', sprnum)
f = HISTORY ('SF', sprnum)

def lin0(t): return (-1, 0, 0)
PRESCRIBE (parnum, linear=lin0)

DEM (1.0, h)

#def lin1(t): return (1, 0, 0)
#PRESCRIBE (parnum, linear=lin1)

#DEM (2.0, h)

#from scipy.interpolate import interp1d

# TODO: compare the spline made of the input spring g(x) with the spline made of (x,f) history

# TODO: loop over many such curves and perform comparisons

print 'Generating spring x-f plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (x, f)
  plt.xlabel ('stroke $(m)$')
  plt.ylabel ('force $(N)$')
  plt.savefig ('tests/spring_x_f.png')

except:
  print 't = ', x
  print 'z = ', f
