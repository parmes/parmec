# PARMEC test --> TORSION_SPRING command test

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

# roll = global z
#TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 0, 1))
# roll = global x
#TORSION_SPRING (parnum, -1, (0, 1, 0), (1, 0, 0), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(1, 0, 0))
# roll = global y
#TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4]) #, droll=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 1, 0))

# yaw = global z
#TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 1, 0), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 0, 1))
# yaw = global x
#TORSION_SPRING (parnum, -1, (0, 1, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(1, 0, 0))
# yaw = global y
#TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4]) #, dyaw=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 1, 0))

# pitch = global z
#TORSION_SPRING (parnum, -1, (0, 0, 1), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 0, 1))
# pitch = global x
#TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(1, 0, 0))
# pitch = global y
#TORSION_SPRING (parnum, -1, (0, 1, 0), (0, 0, 1), kpitch=[-1,1E4, 1,-1E4]) #, dpitch=[-1, 5E2, 1, -5E2])
#VELOCITY (parnum, angular=(0, 1, 0))

# yaw = global z; rotate about all three axes with different damping rates
TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
                                                  droll=[-1, 7E2, 1, -7E2], dyaw=[-1, 5E2, 1, -5E2], dpitch=[-1, 3E2, 1, -3E2])
VELOCITY (parnum, angular=(1, 1, 1))

t = HISTORY ('TIME')
ox = HISTORY ('OX', parnum)
oy = HISTORY ('OY', parnum)
oz = HISTORY ('OZ', parnum)

h = 0.001

print 'Time step:', h

DEM (5.0, h, (0.05, h))

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
  plt.savefig ('tests/torsion_spring_oxyz.png')

except:
  print 't = ', t
  print 'ox = ', ox
  print 'oy = ', oy
  print 'oz = ', oz
