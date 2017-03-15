# LS-DYNA keyfile to PARMEC input file generator
import sys
import time
sys.path.append('.')
from keyfileparse import Keyfile
from sets import Set

# start timer
tic = time.clock()

# check input syntax
if len(sys.argv) < 3:
  print 'SYNOPSIS: python [--skip_general_nonlinear_springs] keyfile2parmec.py path_to_keyfile path_to_parmec_file'
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

# parse input parameters
skip_general_nonlinear_springs = 0
key_path_index = 1
out_path_index = 2
for i in range (1,len(sys.argv)):
  if sys.argv[i] == '--skip_general_nonlinear_springs':
    skip_general_nonlinear_springs = 1
  elif sys.argv[i].endswith('.key'):
    key_path_index = i
  elif sys.argv[i].endswith('.py'):
    out_path_index = i

# begin processing files
print 'Parsing keyfile...'
keyfile = Keyfile(sys.argv[key_path_index])

print 'Writing parmec file (rigid bodies)...'
parmec = open(sys.argv[out_path_index], 'w')

parmec.write ('pid2num = {} # PART_INERTIA to particle number mapping\n')
parmec.write ('eid2num = {} # ELEMENT_DISCRETE to spring number mapping\n')
mcnod2pid = {} # mass center node to part inertia mapping

# define rigid bodies
parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# rigid bodies\n')
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

  mcnode = pi['NODEID']
  center = keyfile.NODES[mcnode]

  parmec.write ('num = ANALYTICAL (' +\
                'inertia = [%g, %g, %g, %g, %g, %g], ' %\
		(I_glo[0], I_glo[4], I_glo[8], I_glo[1], I_glo[2], I_glo[3]) +\
		'mass = %g, ' % pi['TM'] +\
		'position = (%g, %g, %g)' % (center[0], center[1], center[2])+\
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
  mcnod2pid[mcnode] = pi['PID']

parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# prescribe acceleration\n')
parmec.write ('#\n')

#DEFINE_CURVE
print 'Writing parmec file (load curves)...'
parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# spring curves\n')
parmec.write ('#\n')
for lc in keyfile['DEFINE_CURVE']:
  curve = []
  for (t,v) in zip(lc['A1'], lc['O1']):
    curve.append (t)
    curve.append (v)
  parmec.write ('curve%d = %s\n' % (lc['LCID'], str(curve)))
  parmec.write ('tms%d = TSERIES(curve%d)\n' % (lc['LCID'], lc['LCID']))
parmec.write ('curve0 = [-1.0, 0.0, 1.0, 0.0]\n')
parmec.write ('tms0 = TSERIES(0.0)\n\n')

#BOUNDARY_PRESCRIBED_MOTION_NODE
print 'Writing parmec file (boundary conditions)...'
if keyfile.getcard('BOUNDARY_PRESCRIBED_MOTION_NODE') != None:
  pidset = Set()
  parmec.write ('# default zero signals\n')
  for bc in keyfile['BOUNDARY_PRESCRIBED_MOTION_NODE']:
    pid = mcnod2pid[bc['NID']]
    if pid == None:
      print 'ERROR: while prescribing nodal motion -->'
      print '       node with ID =', bc['NID'], 'was not defined as a mass center in a PART_INERTIA card'
      sys.exit(1)
    if pid not in pidset:
      parmec.write ('ACC%d_LIN_X = tms0\n' % pid)
      parmec.write ('ACC%d_LIN_Y = tms0\n' % pid)
      parmec.write ('ACC%d_LIN_Z = tms0\n' % pid)
      parmec.write ('ACC%d_ANG_X = tms0\n' % pid)
      parmec.write ('ACC%d_ANG_Y = tms0\n' % pid)
      parmec.write ('ACC%d_ANG_Z = tms0\n' % pid)
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
    if dof == 1:
      parmec.write ('ACC%d_LIN_X = tms%d\n' % (pid, lcid))
    elif dof == 2:
      parmec.write ('ACC%d_LIN_Y = tms%d\n' % (pid, lcid)) 
    elif dof == 3:
      parmec.write ('ACC%d_LIN_Z = tms%d\n' % (pid, lcid))
    elif dof == 5:
      parmec.write ('ACC%d_ANG_X = tms%d\n' % (pid, lcid))
    elif dof == 6:
      parmec.write ('ACC%d_ANG_Y = tms%d\n' % (pid, lcid))
    elif dof == 7:
      parmec.write ('ACC%d_ANG_Z = tms%d\n' % (pid, lcid))

  parmec.write ('\n')
  parmec.write ('# prescribe acceleration\n')
  for pid in pidset:
    parmec.write ('ACC%d_LIN = (ACC%d_LIN_X, ACC%d_LIN_Y, ACC%d_LIN_Z)\n' % (pid, pid, pid, pid))
    parmec.write ('ACC%d_ANG = (ACC%d_ANG_X, ACC%d_ANG_Y, ACC%d_ANG_Z)\n' % (pid, pid, pid, pid))
    parmec.write ('PRESCRIBE (pid2num[%d], ACC%d_LIN, ACC%d_ANG, "'"aa"'")\n' % (pid, pid, pid))
    parmec.write ('\n')

#BOUNDARY_PRESCRIBED_MOTION_SET
if keyfile.getcard('BOUNDARY_PRESCRIBED_MOTION_SET') != None:
  pidset = Set()
  parmec.write ('# default zero signals\n')
  for bc in keyfile['BOUNDARY_PRESCRIBED_MOTION_SET']:
    for nid in keyfile.SET_NODE_LISTS[bc['SID']]:
      pid = mcnod2pid[nid]
      if pid == None:
	print 'ERROR: while prescribing nodal motion -->'
	print '       node with ID =', nid, 'was not defined as a mass center in a PART_INERTIA card'
	sys.exit(1)
      if pid not in pidset:
	parmec.write ('ACC%d_LIN_X = tms0\n' % pid)
	parmec.write ('ACC%d_LIN_Y = tms0\n' % pid)
	parmec.write ('ACC%d_LIN_Z = tms0\n' % pid)
	parmec.write ('ACC%d_ANG_X = tms0\n' % pid)
	parmec.write ('ACC%d_ANG_Y = tms0\n' % pid)
	parmec.write ('ACC%d_ANG_Z = tms0\n' % pid)
      pidset.add (pid)

  parmec.write ('\n')
  parmec.write ('# perscribed signals\n')
  for bc in keyfile['BOUNDARY_PRESCRIBED_MOTION_SET']:
    dof = bc['DOF']
    vad = bc['VAD']
    lcid = bc['LCID']
    for nid in keyfile.SET_NODE_LISTS[bc['SID']]:
      pid = mcnod2pid[nid]
      if dof not in (1, 2, 3, 5, 6, 7):
	print 'ERROR: prescribed node motion DOF not in (1, 2, 3, 5, 6, 7) set'
	sys.exit(1)
      if vad != 1:
	print 'ERROR: prescribed node motion VAD != 1 (acceleration)'
	sys.exit(1)
      if dof == 1:
	parmec.write ('ACC%d_LIN_X = tms%d\n' % (pid, lcid))
      elif dof == 2:
	parmec.write ('ACC%d_LIN_Y = tms%d\n' % (pid, lcid))
      elif dof == 3:
	parmec.write ('ACC%d_LIN_Z = tms%d\n' % (pid, lcid))
      elif dof == 5:
	parmec.write ('ACC%d_ANG_X = tms%d\n' % (pid, lcid))
      elif dof == 6:
	parmec.write ('ACC%d_ANG_Y = tms%d\n' % (pid, lcid))
      elif dof == 7:
	parmec.write ('ACC%d_ANG_Z = tms%d\n' % (pid, lcid))

    parmec.write ('\n')
    parmec.write ('# prescribe acceleration\n')
    for pid in pidset:
      parmec.write ('ACC%d_LIN = (ACC%d_LIN_X, ACC%d_LIN_Y, ACC%d_LIN_Z)\n' % (pid, pid, pid, pid))
      parmec.write ('ACC%d_ANG = (ACC%d_ANG_X, ACC%d_ANG_Y, ACC%d_ANG_Z)\n' % (pid, pid, pid, pid))
      parmec.write ('PRESCRIBE (pid2num[%d], ACC%d_LIN, ACC%d_ANG, "'"aa"'")\n' % (pid, pid, pid))
      parmec.write ('\n')

print "Mapping keyfile nodes to part inertia cards..."
nod2pid = {} # node to part inertia mapping
for cxns in keyfile['CONSTRAINED_EXTRA_NODES_SET']:
  pid = cxns['PID']
  sid = cxns['NSID']
  for i in keyfile.SET_NODE_LISTS[sid]:
    nod2pid[i] = pid

print "Mapping keyfile parts to material identifiesrs..."
pid2mid = {} # part to material mapping
for part in keyfile['PART']:
  pid = part['PID']
  mid = part['MID']
  pid2mid[pid] = mid

print "Mapping orientation vectors..."
vid2dso = {} # vector id to orientation card mapping
for dso in keyfile['DEFINE_SD_ORIENTATION']:
  vid2dso[dso['VID']] = dso

print "Mapping materials..."
mid2mat = {} # material id to material mapping
for mat in keyfile['MAT_SPRING_NONLINEAR_ELASTIC']:
  mid2mat[mat['MID']] = mat

#ELEMENT_DISCRETE
print 'Writing parmec file (springs)...'
parmec.write ('\n')
parmec.write ('#\n')
parmec.write ('# springs\n')
parmec.write ('#\n')
for ed in keyfile['ELEMENT_DISCRETE']:
  n1 = ed['N1']
  n2 = ed['N2']
  pid1 = nod2pid[n1]
  pid2 = nod2pid[n2]
  if pid1 == pid2:
    print "WARNING: spring nodes (%d, %d) belong to the same part %d;\n" % (n1, n2, pid1)
    print "         -> skipping this ELEMENT_DISCRETE card;\n"
    continue
  if n1 == None or n2 == None:
    print 'ERROR: invalid node to part inertia cards mapping'
    sys.exit(1)
  pnt1 = keyfile.NODES[n1]
  pnt2 = keyfile.NODES[n2]

  vid = ed['VID']
  if vid != 0:
    dso = vid2dso[vid]
    if dso == None:
      print 'ERROR: did not find DEFINE_SD_ORIENTATION card with VID = ', vid
      sys.exit(1)
    if dso['IOP'] == 0:
      planar = 'OFF'
    elif dso['IOP'] == 1:
      planar = 'ON'
    else:
      print 'ERROR: unsupported IOP != (0 or 1) in DEFINE_SD_ORIENTATION card with VID = ', vid
      sys.exit(1)
    direct = (dso['XT'], dso['YT'], dso['ZT']) # fixed direction
  else: direct = None # direction = (pnt2-pnt1)/|pnt2-pnt1|

  pid = ed['PID']
  mid = pid2mid[pid]
  try:
    mat = mid2mat[mid] # this can give rise to exception if 'mid' is not MAT_SPRING_NONLINEAR_ELASTIC
    lcd = mat['LCD']
    lcr = mat['LCR']
    spring = 'curve%d' % lcd
    damper = 'curve%d' % lcr
    if direct != None:
      parmec.write ('num = SPRING (pid2num[%d], %s, pid2num[%d], %s, %s, %s, %s, "'"%s"'")\n' % \
		   (pid1, str(pnt1), pid2, str(pnt2), spring, damper, direct, planar))
    else:
      parmec.write ('num = SPRING (pid2num[%d], %s, pid2num[%d], %s, %s, %s)\n' % \
		   (pid1, str(pnt1), pid2, str(pnt2), spring, damper))
    parmec.write ('eid2num[%d] = num\n' % ed['EID'])
  except:
    if not skip_general_nonlinear_springs:
      mat = keyfile.getcard('MAT_SPRING_GENERAL_NONLINEAR', MID=mid) # TODO: optimize out by mapping
      if mat != None:
	lcdl = mat['LCDL']
	lcdu = mat['LCDU']
	spring = 'curve%d' % lcdl
	unload = 'curve%d' % lcdu
	if direct != None:
	  parmec.write ('num = SPRING (pid2num[%d], %s, pid2num[%d], %s, spring=%s, direction=%s, planar="'"%s"'", unload=%s)\n' % \
		       (pid1, str(pnt1), pid2, str(pnt2), spring, direct, planar, unload))
	else:
	  parmec.write ('num = SPRING (pid2num[%d], %s, pid2num[%d], %s, spring=%s, unload=%s)\n' % \
		       (pid1, str(pnt1), pid2, str(pnt2), spring, unload))
	parmec.write ('eid2num[%d] = num\n' % ed['EID'])
      else:
	print 'ERROR: MAT_SPRING_GENERAL_NONLINEAR with MID', mid, 'has not been found'

lbz = keyfile.getcard('LOAD_BODY_Z')
if lbz != None:
  print 'Writing parmec file (gravity)...'
  parmec.write ('#\n')
  parmec.write ('# gravity\n')
  parmec.write ('#\n')
  lcid = lbz['LCID']
  parmec.write ('gz = TSERIES(zip(curve%d[::2], [-g for g in curve%d[1::2]]))\n' % (lcid, lcid))
  parmec.write ('GRAVITY (0.0, 0.0, gz)\n')

gd = keyfile.getcard('DAMPING_GLOBAL')
if gd != None:
  print 'Writing parmec file (global damping)...'
  parmec.write ('#\n')
  parmec.write ('# global damping\n')
  parmec.write ('#\n')
  lcid = gd['LCID']
  parmec.write ('DMP_LIN_X = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['STX'], lcid))
  parmec.write ('DMP_LIN_Y = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['STY'], lcid))
  parmec.write ('DMP_LIN_Z = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['STZ'], lcid))
  parmec.write ('DMP_ANG_X = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['SRX'], lcid))
  parmec.write ('DMP_ANG_Y = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['SRY'], lcid))
  parmec.write ('DMP_ANG_Z = TSERIES(zip(curve%d[::2], [%g*d for d in curve%d[1::2]]))\n' % (lcid, gd['SRZ'], lcid))
  parmec.write ('dmplin = (DMP_LIN_X, DMP_LIN_Y, DMP_LIN_Z)\n');
  parmec.write ('dmpang = (DMP_ANG_X, DMP_ANG_Y, DMP_ANG_Z)\n');
  parmec.write ('DAMPING (dmplin, dmpang)\n')

print 'Writing parmec file (output intervals)...'
parmec.write ('#\n')
parmec.write ('# output intervals\n')
parmec.write ('#\n')
d3plot = keyfile.getcard('DATABASE_BINARY_D3PLOT')
if d3plot != None:
  lcdt = d3plot['LCDT']
  parmec.write ('dt_files = tms%d\n' % lcdt)
else: parmec.write ('dt_fiels = 0\n')

d3thdt = keyfile.getcard('DATABASE_BINARY_D3THDT')
if d3thdt != None:
  lcdt = d3thdt['LCDT']
  parmec.write ('dt_hist = tms%d\n' % lcdt)
else: parmec.write ('dt_hist = 0\n')

ct = keyfile.getcard('CONTROL_TERMINATION')
ts = keyfile.getcard('CONTROL_TIMESTEP')
if ct != None and ts != None:
  print 'Writing parmec file (commented out simulation control)...'
  parmec.write ('#\n')
  parmec.write ('# run simulation\n')
  parmec.write ('#\n')
  parmec.write ('# DEM(%g, %g, (dt_files, dt_hist))\n' % (ct['ENDTIM'], ts['DTINIT']))

parmec.close()

# end timer
toc = time.clock()
print 'Translated', len(keyfile['PART_INERTIA']), 'rigid bodies and', \
       len(keyfile['ELEMENT_DISCRETE']), 'springs in ', toc-tic, 'seconds'
