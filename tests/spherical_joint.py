# PARMEC test --> spherical joint made of SPRING and TORSION_SPRIN

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes1 = [0, 0, 0,
          1, 0, 0,
	  1, 1, 0,
	  0, 1, 0,
	  0, 0, 1,
	  1, 0, 1,
	  1, 1, 1,
	  0, 1, 1]
nodes2 = [0.25, 0.25, 1,
          0.75, 0.25, 1,
	  0.75, 0.75, 1,
	  0.25, 0.75, 1,
	  0.25, 0.25, 4,
	  0.75, 0.25, 4,
	  0.75, 0.75, 4,
	  0.25, 0.75, 4]
elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

part1 = MESH (nodes1, elements, matnum, colors)
part2 = MESH (nodes2, elements, matnum, colors)

# fix part1 at 0-level
p1 = tuple(nodes1[0:3])
p2 = tuple(nodes1[3:6])
p3 = tuple(nodes1[6:9])
p4 = tuple(nodes1[9:12])
SPRING (part1, p1, -1, p1, spring=[-1, -1E6, 1, 1E6], dashpot=1.0)
SPRING (part1, p2, -1, p2, spring=[-1, -1E6, 1, 1E6], dashpot=1.0)
SPRING (part1, p3, -1, p3, spring=[-1, -1E6, 1, 1E6], dashpot=1.0)
SPRING (part1, p4, -1, p4, spring=[-1, -1E6, 1, 1E6], dashpot=1.0)

# spherical joint at point (0.5,0.5,1)
SPRING (part1, (0.5,0.5,1), part2, (0.5,0.5,1), spring=[-1, -1E6, 1, 1E6], dashpot=1.0)
TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0),
  kroll=[-2, 1E5, -1, 0, 1, 0, 2, -1E5], droll=1.0,
  kyaw = [-1, 1E5, 1, -1E5], dyaw = 1,
  kpitch=[-2, 1E5, -1, 0, 1, 0, 2, -1E5], dpitch=1.0)
VELOCITY (part2, angular=(0.1, 0.1, 0))

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
ox = HISTORY ('OX', part2)
oy = HISTORY ('OY', part2)
oz = HISTORY ('OZ', part2)

h = 0.1 * CRITICAL()

print 'Time step:', h

DEM (15.0, h, (0.05, h))

print 'Generating plots ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, ox, label='ox')
  plt.plot (t, oy, label='oy')
  plt.plot (t, oz, label='oz')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time (s)')
  plt.ylabel ('omega (rad/s)')
  plt.legend (loc = 'upper right')
  plt.title ('particle 2')
  plt.savefig ('tests/spherical_joint_oxyz.png')

except:
  print 't = ', t
  print 'ox = ', ox
  print 'oy = ', oy
  print 'oz = ', oz
