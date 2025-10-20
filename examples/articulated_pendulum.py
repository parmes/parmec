# PARMEC example --> articulated pendulum
# using BALL_JOINT and TORSION_SPRING commands

nele=2 # number of elements
tspr=True # enable torsion springs
ktrq=10 # flex stiffness
dtrq=2 # flex damping ratio
step=.0005 # time step

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

for i in range (0, nele):
  nodes = [i,   0,   0,
	   i+1, 0,   0,
	   i+1, 0.1, 0,
	   i,   0.1, 0,
	   i,   0,   0.1,
	   i+1, 0,   0.1,
	   i+1, 0.1, 0.1,
	   i,   0.1, 0.1]

  elements = [8, 0, 1, 2, 3, 4, 5, 6, 7, matnum]

  colors = [1, 4, 0, 1, 2, 3, 2, 4, 4, 5, 6, 7, 3]

  parnum = MESH (nodes, elements, matnum, colors)

  if i: 
    BALL_JOINT (parnum, (i, 0.05, 0.05), prevpar)

    if tspr:
      TORSION_SPRING (prevpar, parnum, (1, 0, 0), (0, 0, 1),
        # (roll, pitch) cone spring and damper curves restrain
        # freedom of rotationally invariant bending at joints
        kroll=[-1, ktrq, -0.5, 0, 0.5, 0, 1, -ktrq], droll=dtrq,
        kyaw=[-1, ktrq, 1, -ktrq], dyaw=dtrq, # blocked twist
        cone = ('roll', 'pitch'),
        refpnt = (i, 0.05, 0.05))

  prevpar = parnum

BALL_JOINT (parnum, (i+1, 0.05, 0.05))

if tspr:
  TORSION_SPRING (parnum, -1, (1, 0, 0), (0, 0, 1),
    kyaw=[-1, ktrq, 1, -ktrq], dyaw=dtrq, # blocked twist 
    refpnt = (i+1, 0.05, 0.05))

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
px0 = HISTORY ('PX', 0, (1, 0.05, 0.05))
py0 = HISTORY ('PY', 0, (1, 0.05, 0.05))
pz0 = HISTORY ('PZ', 0, (1, 0.05, 0.05))
if nele > 1:
  px1 = HISTORY ('PX', 1, (1, 0.05, 0.05))
  py1 = HISTORY ('PY', 1, (1, 0.05, 0.05))
  pz1 = HISTORY ('PZ', 1, (1, 0.05, 0.05))
if tspr:
  global tqr, tqp, tqy
  tqr = HISTORY ('TRQTOT_R', 0)
  tqp = HISTORY ('TRQTOT_P', 0)
  tqy = HISTORY ('TRQTOT_Y', 0)

DEM (10., step, (0.05, step))

print('Generating plots ...', end=' ')
import sys
sys.stdout.flush()

try:
  import matplotlib.pyplot as plt

  dp = []
  if nele > 1:
    for (x0,y0,z0,x1,y1,z1) in zip(px0,py0,pz0,px1,py1,pz1):
      dp.append (((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)**0.5)
  else:
    for (x0,y0,z0) in zip(px0,py0,pz0):
      dp.append (((x0-1.0)**2+(y0-0.05)**2+(z0-0.05)**2)**0.5)
  plt.clf ()
  plt.plot (t, dp)
  plt.xlim ((0, t[-1]))
  plt.xlabel ('time (s)')
  plt.ylabel ('|p(body 0) - p(body 1)| (m)')
  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
  plt.title ('Joint point p motion difference (bodies: 0,1)')
  plt.savefig ('examples/articulated_pendulum_dp_n%d_t%d.png' % (nele, int(tspr)))

  if tspr:
    trp = []
    for (r,q) in zip(tqr,tqp): trp.append ((r*r+q*q)**0.5)
    plt.clf ()
    plt.plot (t, trp)
    plt.xlim ((0, t[-1]))
    plt.xlabel ('time (s)')
    plt.ylabel ('$\sqrt{t_{roll}^{2}+t_{pitch}^{2}}$ (Nm)')
    plt.title ('Resultant bending torque')
    plt.savefig ('examples/articulated_pendulum_trp.png')
    plt.clf ()
    plt.plot (t, tqy)
    plt.xlim ((0, t[-1]))
    plt.xlabel ('time (s)')
    plt.ylabel ('$t_{yaw}$ (Nm)')
    plt.title ('Twsiting torque')
    plt.savefig ('examples/articulated_pendulum_tyaw.png')

  print('ok.')

except:
  print('failed. (check matplotlib)')
