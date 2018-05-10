# PARMEC test --> SPRING curve correctness tests
import random

h = 0.02 # time step within the [0,2]s time range
itest = 0 # test counter index
ncol = 3 # number of columns in the juxtaposed PNG image
nrow = 4 # number of rows in the juxtaposed PNG image
convert = '' # juxtaposed PNG convertion command
allpassed = True # all tests passed flag

for irow in range(0,nrow):
  for icol in range(0,ncol):

# reset all data
    RESET()

# default analytical particle
    MATERIAL (1, 1, 0.25)
    parnum = ANALYTICAL ()

# generate random spring curve within the [-1, -1, 1, 1] box
    n = random.randint(2,20)
    t = []
    s = []
    t.append(0.)
    for i in range(1,n):
      dt = random.uniform(1.,10.)
      t.append(t[i-1]+dt)
    for i in range(0,n):
      t[i] = -1.0 + 2.0*t[i]/t[-1]
      s.append(random.uniform(-1.,1.))
    ts = list(sum(zip(t,s),()))
    sprnum = SPRING (parnum, (0, 0, 0), -1, (0, 0, 0), ts, direction = (-1, 0, 0))

# go back by 1m (apply -1m/s velocity for 1s)
    def lin0(t): return (-1, 0, 0)
    PRESCRIBE (parnum, linear=lin0)
    DEM (1.0, h)

# start recording stroke-force time history
    x = HISTORY ('STROKE', sprnum)
    f = HISTORY ('SF', sprnum)

# go forward by 2m: span [-1,1] range
    def lin1(t): return (1, 0, 0)
    PRESCRIBE (parnum, linear=lin1)
    DEM (2.0+h, h)

# test numerical discrepancy between input-output curves
    from scipy.interpolate import interp1d
    import numpy as np
    g0 = interp1d(t,s)
    g1 = interp1d(x,f)
    dgok = True
    eps = 0.1
    for w in np.arange(-1.,x[-1]-h,h):
     dg = abs(g0(w)-g1(w))
     if dg > eps*h:
       print '|input(%g)-output(%g)|=|(%g)-(%g)|=%g > 0.1*h (=%g*%g=%g)' % (w, w, g0(w), g1(w), dg, eps, h, eps*h)
       dgok = False
       break
    print 'Test %d: input-output curve difference test: %s' % (itest, 'PASSED' if dgok else 'FAILED')
    if not dgok: allpassed = False

    try:
      import matplotlib
      matplotlib.rcParams.update({'font.size': 14})
      import matplotlib.pyplot as plt
      pngpath = 'tests/spring_curve_test%d.png' % itest
      print 'Test %d: plotting x-f graps to %s' % (itest, pngpath)
      plt.clf ()
      plt.plot (t, s, label='input')
      plt.plot (x, f, label='output', ls='None', marker='.', markevery=2)
      plt.legend (loc = 'upper right')
      plt.xlabel ('stroke $(m)$')
      plt.ylabel ('force $(N)$')
      plt.savefig (pngpath, bbox_inches='tight', pad_inches=0)
      if itest == 0: convert = 'convert '
    except: pass

# iterate to the next test index
    itest = itest + 1

# append convert command syntax for PNG array juxtaposition
    if convert <> '':
      if icol == 0: convert += '\('
      convert += ' ' + pngpath
      if icol == ncol-1: convert += ' +append \)'
  if convert <> '':
    convert += ' -append '
if convert <> '':
  convert += 'tests/spring_curve_tests.png'

# use Image Magic's convert tool to juxtapose PNG files
import subprocess
print 'Running:', convert
process = subprocess.Popen(convert, shell=True)
process.wait()
cleanup = 'rm tests/spring_curve_test{0..%d}.png' % (itest-1)
print 'Running:', cleanup
process = subprocess.Popen(cleanup, shell=True)
process.wait()

# print final status and info
print '===================================='
print 'spring_curve_tests.py STATUS: %s' % ('PASSED' if allpassed else 'FAILED')
print '===================================='
if convert <> '':
  print '....................................................'
  print 'juxtaposed PNG file at: tests/spring_curve_tests.png'
  print '....................................................'
