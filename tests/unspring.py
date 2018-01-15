# PARMEC test --> UNSPRING command test

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

spri = [-1,-1E7, 1,1E7]
dash = [-1, -8E5, 1, 8E5]

for i in range (0, 10):
  nodes = [i, 0, 0,
	   i+1, 0, 0,
	   i+1, 1, 0,
	   i, 1, 0,
	   i, 0, 1,
	   i+1, 0, 1,
	   i+1, 1, 1,
	   i, 1, 1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  parnum = MESH (nodes, elements, matnum, colors)

  if i: 
    SPRING (parnum, (i, 0, 0), prevpar, (i, 0, 0), spri, dash)
    SPRING (parnum, (i, 1, 0), prevpar, (i, 1, 0), spri, dash)
    SPRING (parnum, (i, 1, 1), prevpar, (i, 1, 1), spri, dash)
    SPRING (parnum, (i, 0, 1), prevpar, (i, 0, 1), spri, dash)

  prevpar = parnum

SPRING (parnum, (i+1, 0, 0), -1, (i+1, 0, 0), spri, dash)
SPRING (parnum, (i+1, 1, 0), -1, (i+1, 1, 0), spri, dash)
SPRING (parnum, (i+1, 1, 1), -1, (i+1, 1, 1), spri, dash)
SPRING (parnum, (i+1, 0, 1), -1, (i+1, 0, 1), spri, dash)

UNSPRING ([16, 17, 18, 19], [16, 17, 18, 19], (None, 4E5))

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)
ss16 = HISTORY ('SS', 16)
ss18 = HISTORY ('SS', 18)

h = 0.05 * CRITICAL()

print 'Time step:', h

DEM (2.0, h, (0.05, h))

print 'Generating time-z(center) plot ...'

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, z)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/unspring_z.png')

  plt.clf ()
  plt.plot (t, ss16, label='ss16')
  plt.plot (t, ss18, label='ss18')
  plt.xlim ((0, t[-1]))
  plt.legend ('bottom right')
  plt.xlabel ('time $(s)$')
  plt.ylabel ('spring state')
  plt.savefig ('tests/unspring_ss.png')

except:
  print 't = ', t
  print 'z = ', z
  print 'ss16 =', ss16
  print 'ss18 =', ss18
