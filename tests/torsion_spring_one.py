# PARMEC test --> TORSION_SPRING command test (one particle)

for case in range (0,11):

  RESET()

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

  if case == 0:
    # roll = global z
    TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 0, 1))
    ending = 'z_roll'
  elif case == 1:
    # roll = global x
    TORSION_SPRING (parnum, -1, (0, 1, 0), (1, 0, 0), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(1, 0, 0))
    ending = 'x_roll'
  elif case == 2:
    # roll = global y
    TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 1, 0))
    ending = 'y_roll'
  elif case == 3:
    # yaw = global z
    TORSION_SPRING (parnum, -1, (0, 0, 1), (0, 1, 0), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 0, 1))
    ending = 'z_yaw'
  elif case == 4:
    # yaw = global x
    TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(1, 0, 0))
    ending = 'x_yaw'
  elif case == 5:
    # yaw = global y
    TORSION_SPRING (parnum, -1, (0, 1, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 1, 0))
    ending = 'y_yaw'
  elif case == 6:
    # pitch = global z
    TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 0, 1))
    ending = 'z_pitch'
  elif case == 7:
    # pitch = global x
    TORSION_SPRING (parnum, -1, (0, 0, 1), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(1, 0, 0))
    ending = 'x_pitch'
  elif case == 8:
    # pitch = global y
    TORSION_SPRING (parnum, -1, (0, 0, 1), (1, 0, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (parnum, angular=(0, 1, 0))
    ending = 'y_pitch'
  elif case == 9:
    # rotate about all three axes with different damping curves
    TORSION_SPRING (parnum, -1, (0, 0, 1), (1, 0, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
                                                droll=[-1, 7E2, 1, -7E2], dyaw=[-1, 5E2, 1, -5E2], dpitch=[-1, 3E2, 1, -3E2])
    VELOCITY (parnum, angular=(1, 1, 1))
    ending = 'xyz_all_dcurves'
  elif case == 10:
    # rotate about all three axes with different damping rates
    TORSION_SPRING (parnum, -1, (0, 0, 1), (1, 0, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
						      droll=0.25, dyaw=0.5, dpitch=1.0)
    VELOCITY (parnum, angular=(1, 1, 1))
    ending = 'xyz_all_drates'

  t = HISTORY ('TIME')
  ox = HISTORY ('OX', parnum)
  oy = HISTORY ('OY', parnum)
  oz = HISTORY ('OZ', parnum)

  h = 0.001

  print 'Running %s case ...' % ending

  DEM (5.0, h, (0.05, h))

  print 'Generating %s plots ...' % ending

  try:
    import matplotlib.pyplot as plt

    plt.clf ()
    plt.plot (t, ox, label='ox')
    plt.plot (t, oy, label='oy')
    plt.plot (t, oz, label='oz')
    plt.xlim ((0, t[-1]))
    plt.xlabel ('time (s)')
    plt.ylabel ('angle (rad/s)')
    plt.legend (loc = 'upper right')
    plt.title (ending)
    plt.savefig ('tests/torsion_spring_one_%s.png' % ending)

  except:
    print 't = ', t
    print 'ox = ', ox
    print 'oy = ', oy
    print 'oz = ', oz
