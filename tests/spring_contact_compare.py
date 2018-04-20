# compare spring_contact_* tests results

t0 = HISTORY ('TIME', h5file='tests/spring_contact_plane0rb.h5')
z0 = HISTORY ('DZ', h5file='tests/spring_contact_plane0rb.h5', h5last=True)
#FIXME: this seems to be reading incorrect time history --> why?
#z1 = HISTORY ('PZ', h5file='tests/spring_contact0rb.h5', h5last=True)

try:
  import matplotlib.pyplot as plt

  plt.clf ()
  plt.plot (t0, z0, label='cplane')
  #plt.plot (t0, z1, label='cspring')
  plt.xlim ((0, t0[-1]))
  plt.legend(loc = 'lower right')
  plt.xlabel ('time $(s)$')
  plt.ylabel ('z(center) $(m)$')
  plt.savefig ('tests/spring_contact_compare_z.png')
except:
  print 'ERROR: matplotlib plotting has failed'
