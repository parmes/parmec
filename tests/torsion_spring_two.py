# PARMEC test --> TORSION_SPRING command test (two particle)

for case in range (0,11):

  RESET()

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

  if case == 0:
    # roll = global z
    TORSION_SPRING (part1, part2, (1, 0, 0), (0, 0, 1), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 0, 1))
    ending = 'z_roll'
  elif case == 1:
    # roll = global x
    TORSION_SPRING (part1, part2, (0, 1, 0), (1, 0, 0), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(1, 0, 0))
    ending = 'x_roll'
  elif case == 2:
    # roll = global y
    TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0), kroll=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 1, 0))
    ending = 'y_roll'
  elif case == 3:
    # yaw = global z
    TORSION_SPRING (part1, part2, (0, 0, 1), (0, 1, 0), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 0, 1))
    ending = 'z_yaw'
  elif case == 4:
    # yaw = global x
    TORSION_SPRING (part1, part2, (1, 0, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(1, 0, 0))
    ending = 'x_yaw'
  elif case == 5:
    # yaw = global y
    TORSION_SPRING (part1, part2, (0, 1, 0), (0, 0, 1), kyaw=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 1, 0))
    ending = 'y_yaw'
  elif case == 6:
    # pitch = global z
    TORSION_SPRING (part1, part2, (1, 0, 0), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 0, 1))
    ending = 'z_pitch'
  elif case == 7:
    # pitch = global x
    TORSION_SPRING (part1, part2, (0, 0, 1), (0, 1, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(1, 0, 0))
    ending = 'x_pitch'
  elif case == 8:
    # pitch = global y
    TORSION_SPRING (part1, part2, (0, 0, 1), (1, 0, 0), kpitch=[-1,1E4, 1,-1E4])
    VELOCITY (part1, angular=(0, 1, 0))
    ending = 'y_pitch'
  elif case == 9:
    # rotate about all three axes with different damping curves
    TORSION_SPRING (part1, part2, (0, 0, 1), (1, 0, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
                                                droll=[-1, 7E2, 1, -7E2], dyaw=[-1, 5E2, 1, -5E2], dpitch=[-1, 3E2, 1, -3E2])
    VELOCITY (part1, angular=(1, 1, 1))
    ending = 'xyz_all_dcurves'
  elif case == 10:
    # rotate about all three axes with different damping rates
    TORSION_SPRING (part1, part2, (0, 0, 1), (1, 0, 0), kroll=[-1,1E4, 1,-1E4], kyaw=[-1,1E4, 1,-1E4], kpitch=[-1,1E4, 1,-1E4],
						      droll=0.25, dyaw=0.5, dpitch=1.0)
    VELOCITY (part1, angular=(1, 1, 1))
    ending = 'xyz_all_drates'

  t = HISTORY ('TIME')
  ox1 = HISTORY ('OX', part1)
  oy1 = HISTORY ('OY', part1)
  oz1 = HISTORY ('OZ', part1)
  ox2 = HISTORY ('OX', part2)
  oy2 = HISTORY ('OY', part2)
  oz2 = HISTORY ('OZ', part2)

  h = 0.001

  print 'Running %s case ...' % ending

  DEM (5.0, h, (0.05, h))

  print 'Generating %s plots ...' % ending

  try:
    import matplotlib.pyplot as plt

    plt.clf ()
    plt.plot (t, ox1, label='ox1')
    plt.plot (t, oy1, label='oy1')
    plt.plot (t, oz1, label='oz1')
    plt.plot (t, ox2, label='ox2')
    plt.plot (t, oy2, label='oy2')
    plt.plot (t, oz2, label='oz2')
    plt.xlim ((0, t[-1]))
    plt.xlabel ('time (s)')
    plt.ylabel ('angle (rad/s)')
    plt.legend (loc = 'upper right')
    plt.title (ending)
    plt.savefig ('tests/torsion_spring_two_%s.png' % ending)

  except:
    print 't = ', t
    print 'ox1 = ', ox1
    print 'oy1 = ', oy1
    print 'oz1 = ', oz1
    print 'ox2 = ', ox2
    print 'oy2 = ', oy2
    print 'oz2 = ', oz2
