# PARMEC test --> HISTORY command test

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

VELOCITY (parnum, linear = (1, 0, 0), angular = (0, 0, 1))

t = HISTORY ('TIME')
dx0 = HISTORY ('DX', parnum)
dx1 = HISTORY ('DX', parnum, (1, 1, 1))
dy1 = HISTORY ('DY', parnum, (1, 1, 1))
vx = HISTORY ('VX', parnum)
oz = HISTORY ('OZ', parnum)

DEM (10, 0.1)

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t, dx0)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dx(center) $(m)$')
  plt.savefig ('tests/history_dx0.png')

  plt.clf ()
  plt.plot (t, dx1)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dx(1,1,1) $(m)$')
  plt.savefig ('tests/history_dx1.png')

  plt.clf ()
  plt.plot (t, dy1)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('dy(1,1,1) $(m)$')
  plt.savefig ('tests/history_dy1.png')

  plt.clf ()
  plt.plot (t, vx)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('vx $(m/s)$')
  plt.savefig ('tests/history_vx.png')

  plt.clf ()
  plt.plot (t, oz)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time $(s)$')
  plt.ylabel ('oz $(rad/s)$')
  plt.savefig ('tests/history_oz.png')

except:
  print 'time = ', t
  print 'dx(center) = ', dx0
  print 'dx(1,1,1) = ', dx1
  print 'dy(1,1,1) = ', dy1
  print 'vx = ', vx
  print 'oz = ', oz
