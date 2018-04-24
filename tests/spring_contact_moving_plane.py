# PARMEC test --> SPRING contact moving plane test

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes0 = [-0.5, -0.5, -0.2,
          1.5, -0.5, -0.2,
	  1.5, 1.5, -0.2,
	  -0.5, 1.5, -0.2,
	  -0.5, -0.5, 0,
	  1.5, -0.5, 0,
	  1.5, 1.5, 0,
	  -0.5, 1.5, 0]

nodes1 = [0, 0, 1,
          1, 0, 1,
	  1, 1, 1,
	  0, 1, 1,
	  0, 0, 2,
	  1, 0, 2,
	  1, 1, 2,
	  0, 1, 2]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

part0 = MESH (nodes0, elements, matnum, colors)
part1 = MESH (nodes1, elements, matnum, colors)

spring = [-1, 1E7, 0, 0, 1, 0]
dashpot = [-1, -8E5, 1, 8E5]

SPRING (part1, (0, 0, 1), part0, [(0, 0, 0), (0, 0, 1)], spring, dashpot)
SPRING (part1, (1, 0, 1), part0, [(1, 0, 0), (0, 0, 1)], spring, dashpot)
SPRING (part1, (1, 1, 1), part0, [(1, 1, 0), (0, 0, 1)], spring, dashpot)
SPRING (part1, (0, 1, 1), part0, [(0, 1, 0), (0, 0, 1)], spring, dashpot)

VELOCITY (part1, linear=(0, 0, -1), angular=(0.25, 0.5, 0))

t = HISTORY ('TIME')
z = HISTORY ('PZ', part1)

h = 0.1 * CRITICAL()

print 'Time step:', h

DEM (5.0, h, (0.05, h))

print 'Generating time-z(center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.title ('spring_contact_moving_plane')
  plt.savefig ('tests/spring_contact_moving_plane_z.png')

except:
  print 't = ', t
  print 'z = ', z
