# PARMEC test --> SPRING command test: chain of cubes with critical damping

matnum = MATERIAL (1E3, 1E9, 0.25)

prevpar = -1

dratio0 = 1./2.**0.5 # these values are hand picked so that per-particle
dratio1 = 0.5 # damping ratios are equal to 1.0

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

  if i: SPRING (parnum, (i, i, i), prevpar, (i, i, i), [-1,-1E7, 1,1E7], dratio0)

  prevpar = parnum

sprnum = SPRING (parnum, (i+1, i+1, i+1), -1, (i+1, i+1, i+1), [-1,-1E7, 1,1E7], dratio1)

GRAVITY (0., 0., -10.)

t = HISTORY ('TIME')
z = HISTORY ('PZ', parnum)

hcri = CRITICAL()

print 'Critical time step estimate:', hcri

(hslst, hplst) = CRITICAL(sprnum+1, parnum)

print 'Per-spring time step list:', hslst
print 'Per-particle time step list:', hplst

h = hplst[0][0]

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
  plt.savefig ('tests/spring_chain_critical_z.png')

except:
  print 't = ', t
  print 'z = ', z
