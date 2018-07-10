# PARMEC test --> BALL_JOINT command test

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

for i in range (0, 2):
  nodes = [i, i, i,
	   i+1, i, i,
	   i+1, i+1, i,
	   i, i+1, i,
	   i, i, i+1,
	   i+1, i, i+1,
	   i+1, i+1, i+1,
	   i, i+1, i+1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  parnum = MESH (nodes, elements, matnum, colors)

  if i: BALL_JOINT (parnum, (i, i, i), prevpar)

  prevpar = parnum

sprnum = BALL_JOINT (parnum, (i+1, i+1, i+1))

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)

h = 0.001
DEM (10.0, h, (0.05, h))

print 'Generating time-z(center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/ball_joint_chain_z.png')

except:
  print 't = ', t
  print 'z = ', z
