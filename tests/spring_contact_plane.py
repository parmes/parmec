# PARMEC test --> SPRING contact plane test (compare with spring_contact_plane.py)

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes = [0, 0, 1,
         1, 0, 1,
	 1, 1, 1,
	 0, 1, 1,
	 0, 0, 2,
	 1, 0, 2,
	 1, 1, 2,
	 0, 1, 2]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

parnum = MESH (nodes, elements, matnum, colors)

dashpot = [-1, -8E5, 1, 8E5]
spring = [-1, 1E7, 0, 0, 1, 0] # spring engages when contact penetration begins
                               # spring orientation is aligned with plane normal
                               # hence force for negative penetration is positive
			       # to push the point of the other particle away

SPRING (parnum, (0, 0, 1), -1, [(0, 0, 0), (0, 0, 1)], spring, dashpot, friction = 0.1)
SPRING (parnum, (1, 0, 1), -1, [(1, 0, 0), (0, 0, 1)], spring, dashpot, friction = 0.1)
SPRING (parnum, (1, 1, 1), -1, [(1, 1, 0), (0, 0, 1)], spring, dashpot, friction = 0.1)
sprnum = SPRING (parnum, (0, 1, 1), -1, [(0, 1, 0), (0, 0, 1)], spring, dashpot, friction = 0.1)

VELOCITY (parnum, linear=(1.0, 1.0, 0.0), angular=(0.25, 0.5, 1.0))

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)
f = HISTORY ('F', sprnum)
ff = HISTORY ('FF', sprnum)

h = 0.1 * CRITICAL()

print 'Time step:', h

DEM (5.0, h, (0.05, h))

print 'Generating plots ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.title ('spring_contact_plane - PZ')
  plt.savefig ('tests/spring_contact_plane_z.png')

  plt.clf ()
  plt.semilogy()
  plt.plot (t, f, label = 'contact')
  plt.plot (t, ff, label = 'friction')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('force at (0,1,1) $(N)$')
  plt.legend (loc = 'upper right')
  plt.title ('spring_contact_plane - F, FF')
  plt.savefig ('tests/spring_contact_plane_f_ff.png')

except:
  print 't = ', t
  print 'z = ', z
