t = HISTORY ('TIME', h5file = 'examples/spring_system0rb.h5')
dl = HISTORY('|D|', 0, h5file = 'examples/spring_system0rb.h5')
vl = HISTORY('|V|', 0, h5file = 'examples/spring_system0rb.h5', h5last = True)
length = HISTORY('LENGTH', 0, h5file = 'examples/spring_system0sd.h5', h5last = True)
print 'TIME:'
print t
print '|D|:'
print dl
print '|V|:'
print vl
print '|LENGTH|:'
print length
