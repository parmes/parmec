# PARMEC test --> inactive SPRINGS -> UNSPRING activate test

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

spri0 = [-1,1E7, 0, 0, 1, 0]
spri1 = [-1,-1E7, 1,1E7]
dash = [-1, -8E5, 1, 8E5]

sprmap = []
act0 = []

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
    s0 = SPRING (parnum, (i, 0, 0), prevpar, (i, 0, 0), spri1, dash)
    s1 = SPRING (parnum, (i, 1, 0), prevpar, (i, 1, 0), spri1, dash)
    s2 = SPRING (parnum, (i, 1, 1), prevpar, (i, 1, 1), spri1, dash)
    s3 = SPRING (parnum, (i, 0, 1), prevpar, (i, 0, 1), spri1, dash)
    sprmap.append ([s0, s1, s2, s3])

    for j in range(0,24,3): # inactive base plane contact springs -> activated by UNSPRING
      k = SPRING (parnum, tuple(nodes[j:j+3]), -1, [(0, 0, -5), (0, 0, 1)], spri0, dash, inactive = True)
      act0.append (k)

  prevpar = parnum

SPRING (parnum, (i+1, 0, 0), -1, (i+1, 0, 0), spri1, dash)
SPRING (parnum, (i+1, 1, 0), -1, (i+1, 1, 0), spri1, dash)
SPRING (parnum, (i+1, 1, 1), -1, (i+1, 1, 1), spri1, dash)
SPRING (parnum, (i+1, 0, 1), -1, (i+1, 0, 1), spri1, dash)

UNSPRING (sprmap[4], sprmap[4], (None, 4E5), activate = act0)

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)

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
  plt.savefig ('tests/unspring_activate.png')

except:
  print 't = ', t
  print 'z = ', z
