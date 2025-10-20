import shutil

def inplace_change(filename, pairs):
  # Safely read the input filename using 'with'
  with open(filename) as f:
    s = f.read()

  # Safely write the changed content, if found in the file
  with open(filename, 'w') as f:
    for p in pairs:
      s = s.replace(p[0], p[1])
    f.write(s)

print('Generating dynlb4.h, dynlb8.h, condet4.h, condet8.h ...', end=' ')

shutil.copyfile ('parmec.h', 'parmec4.h')
inplace_change ('parmec4.h', [('REAL', 'float')])
shutil.copyfile ('parmec.h', 'parmec8.h')
inplace_change ('parmec8.h', [('REAL', 'double')])

shutil.copyfile ('objs4/condet_ispc.h', 'condet4.h')
inplace_change ('parmec4.h', [('condet_ispc', 'condet4')])
shutil.copyfile ('objs8/condet_ispc.h', 'condet8.h')
inplace_change ('parmec8.h', [('condet_ispc', 'condet8')])

print('done.')
