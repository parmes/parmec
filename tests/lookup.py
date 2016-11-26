# PARMEC test --> LOOKUP table

x = [0.0, 1.0]
y = [0.0, 1.0]
for i in range (0,4096):

  looknum = LOOKUP (x, y)

  if i % 16 == 0:
    print 'Adding table', looknum, '...'

  x = x + [float(i)]
  y = y + [float(i)]
