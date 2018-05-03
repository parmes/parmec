print '(INFO: run examples/spring_system.py prior to this test)'
t = HISTORY ('TIME', h5file = 'examples/spring_system0rb.h5')
px = HISTORY('PX', 0, h5file = 'examples/spring_system0rb.h5')
py = HISTORY('PY', 0, h5file = 'examples/spring_system0rb.h5')
pz = HISTORY('PZ', 0, h5file = 'examples/spring_system0rb.h5')
dl = HISTORY('|D|', 0, h5file = 'examples/spring_system0rb.h5')
vl0 = HISTORY('|V|', 0, h5file = 'examples/spring_system0rb.h5', h5last = True)
vl1 = HISTORY('|V|', 0, point = (px[0]+0.5,py[0]+0.5,pz[0]+0.5), h5file = 'examples/spring_system0rb.h5', h5last = True)
ol = HISTORY('|O|', 0, h5file = 'examples/spring_system0rb.h5', h5last = True)
length = HISTORY('LENGTH', 0, h5file = 'examples/spring_system0sd.h5', h5last = True)
print 'TIME:'
print t
print '|D|:'
print dl
print '|V0|:'
print vl0
print '|V1|:'
print vl1
print '|O|:'
print ol
print '|LENGTH|:'
print length
