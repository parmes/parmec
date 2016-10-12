# LS-DYNA keyfile to PARMEC input file generator
import sys
sys.path.append('.')
from keyfileparse import Keyfile

# check input syntax
if len(sys.argv) < 3:
  print 'SYNOPSIS: python keyfile2parmec.py path_to_keyfile path_to_parmec_file'
  sys.exit(0)

# auxiliary vector/matrix operations
def cross(a, b):
  return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def matmat(a, b):
  return (a[0]*b[0]+a[3]*b[1]+a[6]*b[2],\
          a[1]*b[0]+a[4]*b[1]+a[7]*b[2],\
          a[2]*b[0]+a[5]*b[1]+a[8]*b[2],\
          a[0]*b[3]+a[3]*b[4]+a[6]*b[5],\
          a[1]*b[3]+a[4]*b[4]+a[7]*b[5],\
          a[2]*b[3]+a[5]*b[4]+a[8]*b[5],\
          a[0]*b[6]+a[3]*b[7]+a[6]*b[8],\
          a[1]*b[6]+a[4]*b[7]+a[7]*b[8],\
          a[2]*b[6]+a[5]*b[7]+a[8]*b[8])

def trans(a):
  return (a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8])

# begin processing files
print 'Parsing keyfile...'
keyfile = Keyfile(sys.argv[1])

print 'Writing parmec file (defining rigid bodies)...'
parmec = open(sys.argv[2], 'w')

parmec.write ('pid2num = {} # part inertia to particle mapping\n')
mcnod2pid = {} # mass center node to part inertia mapping

# define rigid bodies
parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# define rigid bodies\n')
parmec.write ('#\n')
for pi in keyfile['PART_INERTIA']:

  if pi['IRCS'] == 1:
    coord = keyfile.getcard ('DEFINE_COORDINATE_SYSTEM', CID=pi['CID'])
    if coord == None:
      print 'ERROR: DEFINE_COORDINATE_SYSTEM card with CID =', pi['CID'], 'was not found'
      sys.exit(1)
    vec_x = (coord['XL'], coord['YL'], coord['ZL'])
    vec_xy = (coord['XP'], coord['YP'], coord['ZP'])
    vec_z = cross(vec_x, vec_xy)
    vec_y = cross(vec_z, vec_x)
    T_mat = vec_x + vec_y + vec_z
    I_loc = (pi['IXX'], pi['IXY'], pi['IXZ'],\
	     pi['IXY'], pi['IYY'], pi['IYZ'],\
	     pi['IXZ'], pi['IYZ'], pi['IZZ'])
    I_glo = matmat(T_mat, matmat(I_loc, trans(T_mat)))
  else:
    I_glo = (pi['IXX'], pi['IXY'], pi['IXZ'],\
	     pi['IXY'], pi['IYY'], pi['IYZ'],\
	     pi['IXZ'], pi['IYZ'], pi['IZZ'])

  parmec.write ('num = ANALYTICAL (' +\
                'inertia = [%g, %g, %g, %g, %g, %g], ' %\
		(I_glo[0], I_glo[4], I_glo[8], I_glo[1], I_glo[2], I_glo[3]) +\
		'mass = %g, ' % pi['TM'] +\
		'position = (%g, %g, %g)' % (pi['XC'], pi['YC'], pi['ZC'])+\
		')\n')

  mat = keyfile.getcard ('MAT_RIGID', MID=str(int(pi['MID']))) # int() due to '01' used sometimes

  if mat == None:
    print 'ERROR: MAT_RIGID card with MID =', int(pi['MID']), 'was not found'
    sys.exit(1)

  con1 = mat['CON1']
  if con1 == 1:
    lincon = '[1.0, 0.0, 0.0]'
  elif con1 == 2:
    lincon = '[0.0, 1.0, 0.0]'
  elif con1 == 3:
    lincon = '[0.0, 0.0, 1.0]'
  elif con1 == 4:
    lincon = '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0]'
  elif con1 == 5:
    lincon = '[0.0, 1.0, 0.0, 0.0, 0.0, 1.0]'
  elif con1 == 6:
    lincon = '[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]'
  elif con1 == 7:
    lincon = '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]'

  con2 = mat['CON2']
  if con2 == 1:
    angcon = '[1.0, 0.0, 0.0]'
  elif con2 == 2:
    angcon = '[0.0, 1.0, 0.0]'
  elif con2 == 3:
    angcon = '[0.0, 0.0, 1.0]'
  elif con2 == 4:
    angcon = '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0]'
  elif con2 == 5:
    angcon = '[0.0, 1.0, 0.0, 0.0, 0.0, 1.0]'
  elif con2 == 6:
    angcon = '[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]'
  elif con2 == 7:
    angcon = '[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]'

  if con1 > 0 and con2 > 0: parmec.write ('CONSTRAIN (num, ' + lincon + ', ' + angcon + ')\n')
  elif con1 > 0: parmec.write ('CONSTRAIN (num, ' + lincon + ')\n')
  elif con2 > 0: parmec.write ('CONSTRAIN (num, angular = ' + angcon + ')\n')

  parmec.write ('pid2num[%d] = num\n' % pi['PID'])
  mcnod2pid[pi['NODEID']] = pi['PID']

parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# prescribe acceleration\n')
parmec.write ('#\n')
parmec.write ('try:\n')
parmec.write ('  from scipy.interpolate import interp1d\n')
parmec.write ('except:\n')
parmec.write ('  print "'"ERROR: SciPy interp1d failed to load -->"'"\n')
parmec.write ('  print "'"       perhaps SciPy needs to be installed"'"\n')
parmec.write ('\n')

from sets import Set
pidset = Set()
parmec.write ('# default zero signals\n')
for bc in keyfile['BOUNDARY_PRESCRIBED_MOTION_NODE']:
  pid = mcnod2pid[bc['NID']]
  if pid == None:
    print 'ERROR: while prescribing nodal motion -->'
    print '       node with ID =', bc['NID'], 'was not defined as a mass center in a PART_INERTIA card'
    sys.exit(1)
  if pid not in pidset:
    parmec.write ('def ACC%d_LIN_X(t): return 0.0\n' % pid)
    parmec.write ('def ACC%d_LIN_Y(t): return 0.0\n' % pid)
    parmec.write ('def ACC%d_LIN_Z(t): return 0.0\n' % pid)
    parmec.write ('def ACC%d_ANG_X(t): return 0.0\n' % pid)
    parmec.write ('def ACC%d_ANG_Y(t): return 0.0\n' % pid)
    parmec.write ('def ACC%d_ANG_Z(t): return 0.0\n' % pid)
  pidset.add (pid)

parmec.write ('\n')
parmec.write ('# perscribed signals\n')
for bc in keyfile['BOUNDARY_PRESCRIBED_MOTION_NODE']:
  pid = mcnod2pid[bc['NID']]
  dof = bc['DOF']
  vad = bc['VAD']
  lcid = bc['LCID']
  if dof not in (1, 2, 3, 5, 6, 7):
    print 'ERROR: prescribed node motion DOF not in (1, 2, 3, 5, 6, 7) set'
    sys.exit(1)
  if vad != 1:
    print 'ERROR: prescribed node motion VAD != 1 (acceleration)'
    sys.exit(1)
  lc = keyfile.getcard('DEFINE_CURVE', LCID=lcid)
  if lc == None:
    print 'ERROR: DEFINE_CURVE card with ID = ', lcid, 'was not found'
    sys.exit(1)
  tt = lc['A1']
  vv = lc['O1']
  if dof == 1:
    parmec.write ('ACC%d_LIN_X = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))
  elif dof == 2:
    parmec.write ('ACC%d_LIN_Y = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))
  elif dof == 3:
    parmec.write ('ACC%d_LIN_Z = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))
  elif dof == 5:
    parmec.write ('ACC%d_ANG_X = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))
  elif dof == 6:
    parmec.write ('ACC%d_ANG_Y = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))
  elif dof == 7:
    parmec.write ('AC%dC_ANG_Z = interp1d(%s, %s)\n' % (pid, str(tt), str(vv)))

parmec.write ('\n')
parmec.write ('# prescribe acceleration\n')
for pid in pidset:
  parmec.write ('def ACC%d_LIN(t):\n' % pid)
  parmec.write ('  return (float(ACC%d_LIN_X(t)), float(ACC%d_LIN_Y(t)), float(ACC%d_LIN_Z(t)))\n' % (pid, pid, pid))
  parmec.write ('def ACC%d_ANG(t):\n' % pid)
  parmec.write ('  return (float(ACC%d_ANG_X(t)), float(ACC%d_ANG_Y(t)), float(ACC%d_ANG_Z(t)))\n' % (pid, pid, pid))
  parmec.write ('\n')
  parmec.write ('PRESCRIBE (pid2num[%d], ACC%d_LIN, ACC%d_ANG, "'"aa"'")\n' % (pid, pid, pid))

print "Mapping keyfile nodes to part inertia cards..."
nod2pid = {} # node to part inertia mapping
for cxns in keyfile['CONSTRAINED_EXTRA_NODES_SET']:
  pid = cxns['PID']
  sid = cxns['NSID']
  for i in keyfile.SET_NODE_LISTS[sid]:
    nod2pid[i] = pid

print 'Writing parmec file (springs)...'
parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# define springs\n')
parmec.write ('#\n')
for ed in keyfile['ELEMENT_DISCRETE']:
  n1 = ed['N1']
  n2 = ed['N2']
  pid1 = nod2pid[n1]
  pid2 = nod2pid[n2]
  if n1 == None or n2 == None:
    print 'ERROR: invalid node to part inertia cards mapping'
    sys.exit(1)
  pnt1 = keyfile.NODES[n1]
  pnt2 = keyfile.NODES[n2]

  vid = ed['VID']
  dso = keyfile.getcard('DEFINE_SD_ORIENTATION', VID=vid)
  if dso == None:
    print 'ERROR: did not find DEFINE_SD_ORIENTATION card with VID = ', vid
    sys.exit(1)
  #FIXME: undo the below
  #if dso['IOP'] != 0:
    #print 'ERROR: unsupported IOP != 0 in DEFINE_SD_ORIENTATION card with VID = ', vid
    #sys.exit(1)
  direct = (dso['XT'], dso['YT'], dso['ZT'])

  pid = ed['PID']
  part = keyfile.getcard('PART', PID=ed['PID'])
  if part == None:
    print 'ERROR: did not find PART card with PID = ', ed['PID']
    sys.exit(1)
  mid = part['MID']
  mat = keyfile.getcard('MAT_SPRING_NONLINEAR_ELASTIC', MID=part['MID'])
  if mat == None:
    print 'ERROR: did not find MAT_SPRING_NONLINEAR_ELASTIC card with MID = ', mid
    sys.exit(1)
  lcd = mat['LCD']
  lc = keyfile.getcard('DEFINE_CURVE', LCID=lcd)
  if lc == None:
    print 'ERROR: DEFINE_CURVE card with ID = ', lcd, 'was not found'
    sys.exit(1)
  spring = []
  for (t,v) in zip(lc['A1'], lc['O1']):
    spring.append (t)
    spring.append (v)
  lcr = mat['LCR']
  if lcr > 0:
    lc = keyfile.getcard('DEFINE_CURVE', LCID=lcr)
    if lc == None:
      print 'ERROR: DEFINE_CURVE card with ID = ', lcr, 'was not found'
      sys.exit(1)
    damper = []
    for (t,v) in zip(lc['A1'], lc['O1']):
      damper.append (t)
      damper.append (v)
  else: 
    damper = [-1, 0, 1, 0]

  parmec.write ('SPRING (pid2num[%d], %s, pid2num[%d], %s, %s, %s, %s)\n' % \
    (pid1, str(pnt1), pid2, str(pnt2), str(spring), str(damper), str(direct)))

parmec.close()

print 'Translated', len(keyfile['PART_INERTIA']), 'rigid bodies and', \
       len(keyfile['ELEMENT_DISCRETE']), 'springs.'

  #TODO: optimize parmec file size by defining all load curves first
  #      end then only refering to them when used
