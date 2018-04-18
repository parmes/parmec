# PARMEC test --> SPRING command test with stroke offset

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

tmsnum = TSERIES ([0, 0, 5, 0, 10, 1, 15, 1])

SPRING (parnum, (1, 1, 1), -1, (1, 1, 1), [-1,-1E7, 1,1E7], [-1, -8E5, 1, 8E5], offset = tmsnum)

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)

h = 0.3 * CRITICAL()

print 'Time step:', h

DEM (15.0, h, (0.05, h))

print 'Generating time-z(center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/spring_stroke_offset_z.png')

except:
  print 't = ', t
  print 'z = ', z
