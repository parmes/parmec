# PARMEC test --> TORSION_SPRING command test -- two particles

matnum = MATERIAL (1E3, 1E9, 0.25)

nodes1 = [0, 0, 0,
          1, 0, 0,
	  1, 1, 0,
	  0, 1, 0,
	  0, 0, 1,
	  1, 0, 1,
	  1, 1, 1,
	  0, 1, 1]

nodes2 = [0, 0, 1,
          1, 0, 1,
	  1, 1, 1,
	  0, 1, 1,
	  0, 0, 2,
	  1, 0, 2,
	  1, 1, 2,
	  0, 1, 2]

elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

part1 = MESH (nodes1, elements, matnum, colors)

part2 = MESH (nodes2, elements, matnum, colors)

# roll = global z
#TORSION_SPRING (part1, part2, (1, 0, 0), (0, 0, 1), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(0, 0, 1))
# roll = global x
#TORSION_SPRING (part1, part2, (0, 1, 0), (1, 0, 0), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(1, 0, 0))
# roll = global y
#TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(0, 1, 0))

# yaw = global z
TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
VELOCITY (part1, angular=(0, 0, 1))
# yaw = global x
#TORSION_SPRING (part1, part2, (0, 1, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(1, 0, 0))
# yaw = global y
#TORSION_SPRING (part1, part2, (1, 0, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(0, 1, 0))

# pitch = global z
#TORSION_SPRING (part1, part2, (0, 0, 1), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(0, 0, 1))
# pitch = global x
#TORSION_SPRING (part1, part2, (1, 0, 0), (0, 0, 1), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(1, 0, 0))
# pitch = global y
#TORSION_SPRING (part1, part2, (0, 1, 0), (0, 0, 1), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (part1, angular=(0, 1, 0))

# yaw = global z; rotate about all three axes with different damping rates
#TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
#                                                    droll=[-1, 7E2, 1, -7E2], dyaw=[-1, 5E2, 1, -5E2], dpitch=[-1, 3E2, 1, -3E2])
#VELOCITY (part1, angular=(1, 1, 1))

t = HISTORY ('TIME')
ox1 = HISTORY ('OX', part1)
oy1 = HISTORY ('OY', part1)
oz1 = HISTORY ('OZ', part1)
ox2 = HISTORY ('OX', part2)
oy2 = HISTORY ('OY', part2)
oz2 = HISTORY ('OZ', part2)

h = 0.001

print 'Time step:', h

DEM (5.0, h, (0.05, h))

print 'Generating plots ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, ox1, label='ox')
  plt.plot (t, oy1, label='oy')
  plt.plot (t, oz1, label='oz')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time (s)')
  plt.ylabel ('omega (rad/s)')
  plt.legend (loc = 'upper right')
  plt.title ('particle 1')
  plt.savefig ('tests/torsion_spring_oxyz1.png')

  plt.clf ()
  plt.plot (t, ox2, label='ox')
  plt.plot (t, oy2, label='oy')
  plt.plot (t, oz2, label='oz')
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time (s)')
  plt.ylabel ('omega (rad/s)')
  plt.legend (loc = 'upper right')
  plt.title ('particle 2')
  plt.savefig ('tests/torsion_spring_oxyz2.png')

except:
  print 't = ', t
  print 'ox1 = ', ox1
  print 'oy1 = ', oy1
  print 'oz1 = ', oz1
  print 'ox2 = ', ox2
  print 'oy2 = ', oy2
  print 'oz2 = ', oz2
