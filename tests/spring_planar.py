# PARMEC test --> SPRING command test

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

RESTRAIN (parnum, [0, 0, 1], [1, 0, 0, 0, 1, 0, 0, 0, 1])

SPRING (parnum, (.5, .5, 0), -1, (.5, .5, 0), [-1,-1E5, 1,1E5], [-1, -8E3, 1, 8E3], (0, 0, 1), 'ON')

VELOCITY (parnum, (.5, 1, 0))

t = HISTORY ('TIME')
dx = HISTORY ('DX', parnum)
dy = HISTORY ('DY', parnum)

h = 0.3 * CRITICAL()

print 'Time step:', h

DEM (5.0, h, (0.05, h))

print 'Generating time-(dx,dy) (center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, dx, label='dx')
  plt.plot (t, dy, label='dx')
  plt.legend (loc = 'upper right')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/spring_planar_dxdy.png')

except:
  print 't = ', t
  print 'z = ', z
