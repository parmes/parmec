# PARMEC test --> SPRING command test

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

for i in range (0, 10):
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

  if i: SPRING (parnum, (i, i, i), prevpar, (i, i, i), [-1,-1E7, 1,1E7], [-1, -8E5, 1, 8E5])

  prevpar = parnum

sprnum = SPRING (parnum, (i+1, i+1, i+1), -1, (i+1, i+1, i+1), [-1,-1E7, 1,1E7], [-1, -8E5, 1, 8E5])

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)

hcri = CRITICAL()

print 'Critical time step estimate:', hcri

hlst = CRITICAL(sprnum+1)

print 'Per-spring time step list:', hlst

h = 0.1*hcri

print 'Adopted time step:', h

DEM (10.0, h, (0.05, h))

print 'Generating time-z(center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/spring_chain_z.png')

except:
  print 't = ', t
  print 'z = ', z
