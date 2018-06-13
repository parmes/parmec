# PARMEC test --> spherical joint made of SPRING and TORSION_SPRING

# bulk material
matnum = MATERIAL (1E3, 1E9, 0.25)

# two particles
nodes1 = [0, 0, 0,
          1, 0, 0,
	  1, 1, 0,
	  0, 1, 0,
	  0, 0, 1,
	  1, 0, 1,
	  1, 1, 1,
	  0, 1, 1]
nodes2 = [0.25, 0.25, 1,
          0.75, 0.25, 1,
	  0.75, 0.75, 1,
	  0.25, 0.75, 1,
	  0.25, 0.25, 3,
	  0.75, 0.25, 3,
	  0.75, 0.75, 3,
	  0.25, 0.75, 3]
elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]
colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]
part1 = MESH (nodes1, elements, matnum, colors)
part2 = MESH (nodes2, elements, matnum, colors)

# prescribe part1 motion
v0 = 1.
vx = TSERIES ([0, v0, 2.49, v0, 2.50, -v0, 4.99, -v0, 5.0,  v0, 7.49, v0, 7.5, -v0, 9.00, -v0, 10, 0, 20, 0])
vy = TSERIES ([0, 0, 9.99, 0, 10, v0, 12.49, v0, 12.5, -v0, 14.99, -v0, 15, v0, 17.49, v0, 17.5, -v0, 20, -v0])
vz = TSERIES ([0, 0, 20, 0])
PRESCRIBE (part1, linear=(vx, vy, vz))
RESTRAIN (part1, angular=[1, 0, 0, 0, 1, 0, 0, 0, 1])

# spherical joint at point (0.5,0.5,1) co-rotating with part2
klin=5E6 # linear stiffness
dlin=1 # linear damping ratio
ktrq=5E6 # angular stiffness
dtrq=2 # angular damping ratio
SPRING (part1, (0.5,0.5,1), part2, (0.5,0.5,1), spring=[-1, -klin, 1, klin], dashpot=dlin)
TORSION_SPRING (part2, part1, (0, 0, 1), (1, 0, 0),
  kroll=[-1, ktrq, -0.25, 0, 0.25, 0, 1, -ktrq], droll=dtrq, # allow for some freedom ... (--> cone)
  kyaw=[-1, ktrq, 1, -ktrq], dyaw=dtrq, # block yaw rotation
  cone = ('roll', 'pitch')) # ... in the (roll, pitch) space

# apply gravity
GRAVITY (0., 0., -10.)

# save angular velocity histories on the go
t = HISTORY ('TIME')
ox = HISTORY ('OX', part2)
oy = HISTORY ('OY', part2)
oz = HISTORY ('OZ', part2)

# run at half critical step
h = 0.5 * CRITICAL()
print 'Time step:', h

# simulate 20s, outputting XDMF files at 10fps and sample history at 1000fps
DEM (20., h, (0.1, 0.001))

# plot angular velocities if possible
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
  plt.title ('particle 2')
  plt.savefig ('tests/spherical_joint_oxyz.png')
except:
  print 't = ', t
  print 'ox = ', ox
  print 'oy = ', oy
  print 'oz = ', oz
